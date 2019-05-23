import itertools
import networkx as nx
from collections import defaultdict
import numpy as np 
import pandas as pd
import collections 

def major_minor_allele(counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3] }
    l = sorted(d, key=d.get, reverse=True)
    return l[0], l[1]


def build_SNP_linkage_network(read2snps):
    '''
    Calculates the SNP linkage network
    Arguments:
        read2snps = read -> SNPs
    saves it as a dictionary of edges in self.snv_net.
    Writes it out to a file, genome.net for reading in through other programs like node2vec
    '''
    G=nx.Graph()
    G_pos = nx.Graph()
    for read in read2snps:
        snps = read2snps[read]
        for pair in itertools.combinations(snps, 2):
            if G.has_edge(pair[0], pair[1]):
                G[pair[0]][pair[1]]['weight'] += 1
            else:
                G.add_edge(pair[0], pair[1], weight = 1)
                if not G_pos.has_edge(pair[0][0], pair[1][0]):
                    G_pos.add_edge(pair[0][0], pair[1][0])    
    return G, G_pos # haplotype linkage, site_linkage_network

def calculate_ld(snp_linkage_network, site_linkage_network, scaffold, snp_counts, min_ld):
    '''
    Calculates Linkage Disequilibrium for all SNPs in a window.
    '''
    r2_table = defaultdict(list)

    # This just gives potisions
    for edge in site_linkage_network.edges():
        pos_a = edge[0]
        pos_b = edge[1]

        distance = abs(pos_a - pos_b)
        allele_A, allele_a = major_minor_allele(snp_counts[pos_a])
        allele_B, allele_b = major_minor_allele(snp_counts[pos_b])
        # Get frequencies of linkages
        countAB, countAb, countaB, countab  = 0,0,0,0
        if snp_linkage_network.has_edge(  (pos_a, allele_A), (pos_b, allele_B)):
            countAB = snp_linkage_network[(pos_a, allele_A)][(pos_b, allele_B)]['weight']
        if snp_linkage_network.has_edge(  (pos_a, allele_A), (pos_b, allele_b)):
            countAb = snp_linkage_network[(pos_a, allele_A)][(pos_b, allele_b)]['weight']
        if snp_linkage_network.has_edge(  (pos_a, allele_a), (pos_b, allele_B)):
            countaB = snp_linkage_network[(pos_a, allele_a)][(pos_b, allele_B)]['weight']
        if snp_linkage_network.has_edge(  (pos_a, allele_a), (pos_b, allele_b)):
            countab = snp_linkage_network[(pos_a, allele_a)][(pos_b, allele_b)]['weight']
        ld_result = _calc_ld_single(countAB, countAb, countaB, countab, min_ld)
        if ld_result:
            r2_table['scaffold'].append(scaffold)
            r2_table['position_A'].append(pos_a)
            r2_table['position_B'].append(pos_b)
            r2_table['distance'].append(distance)

            for att in ld_result.keys():
                r2_table[att].append(ld_result[att])

            r2_table['A_allele'].append(allele_A)
            r2_table['a_allele'].append(allele_a)
            r2_table['B_allele'].append(allele_B)
            r2_table['b_allele'].append(allele_b)


    # Create r2 linkage table
    r2linkage_table = pd.DataFrame(r2_table)
    return r2linkage_table

def _calc_ld_single(countAB, countAb, countaB, countab, min_ld):
    '''
    A function that calculates the LD between two sites from their genotype counts.
    '''    

    total = countAB + countAb + countaB + countab

    #Requires at least min_ld
    if total > min_ld:

        #calculate allele frequencies
        freq_AB = float(countAB) / total
        freq_Ab = float(countAb) / total
        freq_aB = float(countaB) / total
        freq_ab = float(countab) / total

        freq_A = freq_AB + freq_Ab
        freq_a = freq_ab + freq_aB
        freq_B = freq_AB + freq_aB
        freq_b = freq_ab + freq_Ab


        linkD = freq_AB - freq_A * freq_B

        if freq_a == 0 or freq_A == 0 or freq_B == 0 or freq_b == 0:
            r2 = np.nan
        else:
            r2 = linkD*linkD / (freq_A * freq_a * freq_B * freq_b)

        linkd = freq_ab - freq_a * freq_b

        # calc D-prime
        d_prime = np.nan
        if (linkd < 0):
            denom = max([(-freq_A*freq_B),(-freq_a*freq_b)])
            d_prime = linkd / denom
            
        elif (linkD > 0):
            denom = min([(freq_A*freq_b), (freq_a*freq_B)])
            d_prime = linkd / denom

        ################
        # calc rarefied

        rareify = np.random.choice(['AB','Ab','aB','ab'], replace=True, p=[freq_AB,freq_Ab,freq_aB,freq_ab], size=min_ld)
        freq_AB = float(collections.Counter(rareify)['AB']) / min_ld
        freq_Ab = float(collections.Counter(rareify)['Ab']) / min_ld
        freq_aB = float(collections.Counter(rareify)['aB']) / min_ld
        freq_ab = float(collections.Counter(rareify)['ab']) / min_ld

        freq_A = freq_AB + freq_Ab
        freq_a = freq_ab + freq_aB
        freq_B = freq_AB + freq_aB
        freq_b = freq_ab + freq_Ab

        linkd_norm = freq_ab - freq_a * freq_b

        if freq_a == 0 or freq_A == 0 or freq_B == 0 or freq_b == 0:
            r2_normalized = np.nan
        else:
            r2_normalized = linkd_norm*linkd_norm / (freq_A * freq_a * freq_B * freq_b)

    
        # calc D-prime
        d_prime_normalized = np.nan
        if (linkd_norm < 0):
            denom = max([(-freq_A*freq_B),(-freq_a*freq_b)])
            d_prime_normalized = linkd_norm / denom
            
        elif (linkd_norm > 0):
            denom = min([(freq_A*freq_b), (freq_a*freq_B)])
            d_prime_normalized = linkd_norm / denom

        rt_dict = {}
        for att in ['r2', 'd_prime', 'r2_normalized', 'd_prime_normalized', 'total', 'countAB', \
                    'countAb', 'countaB', 'countab']:
            rt_dict[att] = eval(att)


        return rt_dict

    else:
        return False