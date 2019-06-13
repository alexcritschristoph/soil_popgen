import numpy as np 
import pandas as pd 
from tqdm import tqdm
import allel
import argparse
import os
from collections import defaultdict
import copy
from Bio import SeqIO

from inStrain_lite import generate_snp_model
from inStrain_lite import call_snv_site
from inStrain_lite import SNPprofile
from inStrain_lite.utilities import major_minor_allele

# def calc_fst(allele_freq1, allele_freq2):
#     Hs = allele_freq1*(1-allele_freq1) + allele_freq2 * (1-allele_freq2)
#     avg_freq = np.mean([allele_freq1, allele_freq2])
#     Ht = 2*avg_freq*(1-avg_freq)
#     Fst = (Ht - Hs) / Ht
#     return(Fst)

def create_gene_index(gene_fasta):
    gene_index = []
    complete_genes = 0.0
    partial_genes = 0.0
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)
        if 'partial=00' in record.description:
            complete_genes += 1
            gene_scaf = "_".join(gene.split("_")[:-1])
            # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE 
            gene_start = int(record.description.split("#")[1].strip())-1
            gene_end = int(record.description.split("#")[2].strip())-1
            gene_index.append({"name":gene, 'scaf': gene_scaf, "start": gene_start, "end": gene_end})

        else:
            partial_genes += 1

    print("Notice: " + str(round(partial_genes *100 / (complete_genes + partial_genes))) + "% of genes were incomplete and snps in these genes were marked I. Pi is not calculated for incomplete genes.")
    return gene_index


def main(args):

    ## Step 0: get null model for SNP calling
    null_loc = os.path.dirname(__file__) + '/helper_files/combined_null1000000.txt'
    null_model = generate_snp_model(null_loc)
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}


    ## Step 1: build new counts table from all objects
    s_final = SNPprofile()
    s_final.filename = args.output
    i = 0
    counts_per_block = {}
    s1 = SNPprofile()
    print("loading " + args.input[0])
    s1.load(args.input[0])

    s_final.scaffold_list = s1.scaffold_list
    s_final.counts_table = copy.deepcopy(s1.counts_table)
    
    s2 = SNPprofile()
    print("loading " + args.input[1])
    s2.load(args.input[1])

    for scaf in s2.scaffold_list:
        if scaf not in s_final.scaffold_list:
            sys.exit("Error: scaffold " + scaf + " in " + fn + " not found in initial file. Your inStrain objects were probably not run on the same FASTA.")

    scaf_counter = 0
    for scaf in s2.counts_table:
        s_final.counts_table[scaf_counter] += scaf
        scaf_counter += 1
    i += 1


    # Step 2: call all SNPs for new object
    allele_counts_total  = {}
    allele_counts1 = {}
    allele_counts2 = {}
    snp_table = defaultdict(list)
    scaf_counter = 0
    
    for scaf in tqdm(s_final.counts_table, desc='Calling new SNVs...'):
        pos_counter = 0
        for counts in scaf:
            snp = call_snv_site(counts, min_cov = 5, min_freq = 0.05, model=null_model)

            if snp: # means that there was coverage at this position
                if snp != -1: # means this is a SNP
                    # calculate varBase
                    snp, varbase = major_minor_allele(counts)

                    snp_table['scaffold'].append(s_final.scaffold_list[scaf_counter])
                    snp_table['position'].append(pos_counter)
                    snp_table['varBase'].append(snp)
                    snp_table['conBase'].append(varbase)
                    allele_counts_total[s_final.scaffold_list[scaf_counter] + ":" + str(pos_counter)] = (s_final.counts_table[scaf_counter][pos_counter])
                    allele_counts1[s_final.scaffold_list[scaf_counter] + ":" + str(pos_counter)] = (s1.counts_table[scaf_counter][pos_counter])
                    allele_counts2[s_final.scaffold_list[scaf_counter] + ":" + str(pos_counter)] = (s2.counts_table[scaf_counter][pos_counter])
            pos_counter += 1 # 0 based positions!!
        scaf_counter += 1

    # Step 3: Save new FST_SNP table to disk.
    SNPTable = pd.DataFrame(snp_table)

    FstTable = defaultdict(list)
    for gene in tqdm(create_gene_index(args.gene_file), desc="calculating fst"):
        snps = SNPTable[(SNPTable.scaffold == gene['scaf']) & (SNPTable.position >= gene['start']) & (SNPTable.position <= gene['end'])]
        snp_list = []
        for index, row in snps.iterrows():
            snp_list.append(row['scaffold'] + ":" + str(row['position']))


        # only continue if there are at least 3 snps in this gene
        if len(snp_list) >= 3:
            allele_counts_1 = []
            allele_counts_2 = []
            for snp in snp_list:
                allele_counts_1.append(allele_counts1[snp])
                allele_counts_2.append(allele_counts2[snp])

            allel1 = allel.AlleleCountsArray(allele_counts_1)
            allel2 = allel.AlleleCountsArray(allele_counts_2)
            fst_h = allel.moving_hudson_fst(allel1, allel2, size=len(snp_list))[0] #allel.moving_hudson_fst(a1,a2, size=3)
            nd_1 = np.sum(allel.mean_pairwise_difference(allel1)) / (1 + gene['end'] - gene['start'])
            nd_2 = np.sum(allel.mean_pairwise_difference(allel2)) / (1 + gene['end'] - gene['start'])

            FstTable['gene'].append(gene['name'])
            FstTable['snp_num'].append(len(snp_list))
            FstTable['fst'].append(fst_h)
            FstTable['pi_1'].append(nd_1)
            FstTable['pi_2'].append(nd_2)
            FstTable['cov_1'].append(np.mean(np.sum(allele_counts_1, axis=1)))
            FstTable['cov_2'].append(np.mean(np.sum(allele_counts_2, axis=1)))

    FstTable = pd.DataFrame(FstTable)
    print(np.mean(FstTable['fst']))
    FstTable.to_csv(args.output + '.Fst.tsv', index=False, sep='\t')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """Calculates Fst, nucleotide diversity, tajima's D, Dxy, and delta tajima's D for two populations.""", formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument('input', nargs='+', help="two inStrain object prefixes, space separated. (don't include `.`, only prefixes - must be in the local directory).")
    parser.add_argument("-g", "--gene_file", action="store", required=True, \
        help='Prodigal .fna gene file for this genome ')

    parser.add_argument("-o", "--output", action="store", required=True, \
        help='Output prefix ')

    args = parser.parse_args()
    main(args)
