#!/usr/bin/env python
# inStrain
# from ._version import __version__
import os
import csv
import sys
import pysam
import pickle
import logging
import argparse
import traceback
import itertools
import collections
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from Bio import SeqIO
import concurrent.futures
from subprocess import call
from itertools import permutations
from collections import defaultdict
## local imports
from inStrain_lite.utilities import *
from inStrain_lite.filter_reads import filter_paired_reads
from inStrain_lite.linkage import build_SNP_linkage_network
from inStrain_lite.linkage import calculate_ld

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
C2P = {0:'A', 1:'C', 2:'T', 3:'G'}

###############################################################
######## START MAIN SNP PROFILE CLASS OBJECT ##################
###############################################################

class SNPprofile():
    def __init__(self, **kwargs):
        # parameters
        output = None
        fasta_db = None
        coverage_table = None
        snp_table = None
        linkage_table = None

        # snp calling model 
        null_loc = os.path.dirname(__file__) + '/helper_files/combined_null1000000.txt'
        global null_model
        null_model = generate_snp_model(null_loc)

    def save(self):
        self.coverage_table.to_csv(self.output + '.scaffold_info.tsv', index=False, sep='\t')
        self.snp_table.to_csv(self.output + '.SNVs.tsv', index=False, sep='\t')
        self.linkage_table.to_csv(self.output + '.linkage.tsv', index=False, sep='\t')
        f = open(self.output + ".data", 'wb')
        pickle.dump(self.__dict__, f, 2)
        f.close()

    def load(self, name):
        self.__dict__.clear()
        f = open(name + ".data", 'rb')
        tmp_dict = pickle.load(f)
        f.close()

    def load_fasta(self, args):
        '''
        Load the FASTA sequences to be profiled.
        Return a table listing scaffold name, start, end
        '''
        # PROFILE ALL SCAFFOLDS IN THE .FASTA FILE
        scaff2sequence = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta")) # set up .fasta file
        s2l = {s:len(scaff2sequence[s]) for s in list(scaff2sequence.keys())} # Get scaffold2length
        Fdb = pd.DataFrame(list(s2l.items()), columns=['scaffold', 'end'])
        Fdb['start'] = 0
        self.fasta_db = Fdb
        return True

    def profile(self, bam, filtered_reads, min_freq=0.05, min_ani=0.97, min_mapq=2, min_insert=50, max_insert=1500, min_ld = 30, min_breadth_cov=[0,0], threads=6):
        ''' Profiles one or multiple scaffolds'''
        args = {"min_freq": min_freq, "min_ani": min_ani, "min_mapq":min_mapq, "min_insert":min_insert, "max_insert":max_insert, 'min_ld': min_ld}
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            Sprofiles = [s for s in tqdm(executor.map(self.profile_scaffold,
                                            self.iterate_commands(self.fasta_db, bam, filtered_reads, args)),
                                            desc='Profiling scaffolds: ',
                                            total=len(self.fasta_db['scaffold'].unique()))]
        
        logging.info("Finalizing and saving")

        counts_table = list()
        counts_array = np.zeros(shape=(0, 4), dtype=int)
        coverage_breadth_table = defaultdict(list)
        i = 0
        for Sprof in Sprofiles:

            # create counts table
            counts_table.append(Sprof[0])
            counts_array = np.append(counts_array, Sprof[0], axis=0)
            coverage_breadth_table['scaffold'].append(Sprof[3])
            coverage_breadth_table['coverage'].append(np.mean(np.sum(Sprof[0], axis=1)))
            coverage_breadth_table['breadth'].append((np.sum( Sprof[0], axis=1) > 0).sum() / float( Sprof[0].shape[0]))

            # create snp / linkage tables
            if i == 0:
                snp_table = Sprof[1]
                linkage_table = Sprof[2]
            else:
                snp_table = pd.concat([snp_table, Sprof[1]])
                linkage_table = pd.concat([linkage_table, Sprof[2]])
            i += 1   

        coverage_breadth_table = pd.DataFrame(coverage_breadth_table)
        total_coverage = np.mean(np.sum(counts_array, axis=1))
        total_breadth = (np.sum(counts_array, axis=1) > 0).sum() / float(counts_array.shape[0])
        print(min_breadth_cov)
        total_min_cov_breadth = (np.sum(counts_array, axis=1) >= min_breadth_cov[1]).sum() / float(counts_array.shape[0])

        logging.info("Mean coverage\t" + str(total_coverage))
        logging.info("Mean breadth\t" + str(total_breadth))
        logging.info("Breadth above minimum coverage\t" + str(total_min_cov_breadth))

        if total_min_cov_breadth >= min_breadth_cov[0]:
            logging.info("Passes breadth-cov filter. Will output tables")
            self.counts_table = counts_table
            self.coverage_table = coverage_breadth_table
            self.snp_table = snp_table
            self.linkage_table = linkage_table
            self.save()
        else:
            logging.info("NOTE: does NOT pass breadth-cov filter. Will not output tables")

    def iterate_commands(self, Fdb, bam, filtered_reads, args):
        '''
        Make and iterate profiling commands
        Doing it in this way makes it use way less RAM
        '''
        for index, row in Fdb.iterrows():
            cmd = {"scaffold": row['scaffold'], 'scaf_length': row['end'], "bam": bam, "filtered_reads": filtered_reads, "args":args}
            yield cmd

    def profile_scaffold(self, args):
        snp_table = defaultdict(list)
        read2snps = defaultdict(list)
        samfile = pysam.AlignmentFile(args['bam'])
        try:
            iter = samfile.pileup(args['scaffold'], truncate = True, max_depth=100000,
                                    stepper = 'nofilter', compute_baq= True,
                                    ignore_orphans = True, ignore_overlaps = True,
                                    min_base_quality = 30)
        except ValueError:
            logging.error("ERROR " + args['scaffold'] + " is not in the BAM file.")
            sys.exit(1)

        pileup_counts = np.zeros(shape=(args['scaf_length'], 4), dtype=int) # Holds all pileup counts - alexcc 5/9/2019
        for pileupcolumn in iter:
            counts = get_base_counts(pileupcolumn, args['filtered_reads'])

            counts_sum = np.sum(counts)

            pileup_counts[pileupcolumn.pos,] = counts
            snp = call_snv_site(counts, min_freq = args['args']['min_freq'])
            if snp: # means this is a SNP
                snp, varbase = major_minor_allele(counts)
                # Add to SNP table
                snp_table['scaffold'].append(args['scaffold'])
                snp_table['position'].append(pileupcolumn.pos)
                for b, c in zip(['A', 'C', 'T', 'G'], counts):
                    snp_table[b].append(c)
                snp_table['conBase'].append(snp)
                snp_table['varBase'].append(varbase)
                snp_table['baseCoverage'].append(np.sum(counts))
                # Add to read2snps.
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_name in args['filtered_reads']:
                        read2snps[pileupread.alignment.query_name].append( (pileupcolumn.pos, pileupread.alignment.query_sequence[pileupread.query_position]) )
        samfile.close()

        snp_table = pd.DataFrame(snp_table)
        snp_net, linkage_net = build_SNP_linkage_network(read2snps)
        linkage_table = calculate_ld(snp_net, linkage_net, args['scaffold'], pileup_counts, args['args']['min_ld'])
        return pileup_counts, snp_table, linkage_table, args['scaffold']
 

###############################################################
######## SCRIPT CONTROLLER CLASS  #############################
###############################################################

class Controller():
    '''
    Main controller of the program
    '''

    def main(self, sys_args):
        '''The main method when run on the command line'''
    
        strains = SNPprofile()

        args = parse_arguments(sys_args)
        setup_logger(args.output + '.log')
        args.min_breadth_cov = [int(i) for i in args.min_breadth_cov.split(",")[:2]]
        if not args.output:
            strains.output = args.fasta.split("/")[-1].split(".")[0]
        else:
            strains.output = args.output
        
        strains.load_fasta(args)

        filtered_reads = filter_paired_reads(args.bam, args.fasta, args.min_ani, args.min_mapq, args.min_insert, args.max_insert)

        strains.profile(args.bam, filtered_reads, args.min_freq, args.min_ani, args.min_mapq, args.min_insert, args.max_insert, args.min_ld, args.min_breadth_cov, args.threads)


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="inStrain",
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required positional arguments
    parser.add_argument("bam", help="Sorted .bam file")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")

    # Optional arguments
    parser.add_argument("-o", "--output", action="store", default=None, \
        help='Output prefix')
    parser.add_argument("-p", "--threads", action="store", default=6, type=int, \
        help='Threads to use for multiprocessing')
    parser.add_argument("-s", "--min_ld", action="store", default=20, \
        help='Minimum number of reads connecting two segregating sites to calculate LD between them. (>=)')
    parser.add_argument("-f", "--min_freq", action="store", default=0.05, \
        help='Minimum SNP frequency to confirm a SNP (both this AND the 10^-6 FDR SNP count cutoff must be true).')

    # Read filtering cutoffs
    parser.add_argument("-c", "--min_ani", action="store", default=0.95, type=float, \
        help='Minimum required percent identity of read pairs to consensus. (>=)')
    parser.add_argument("--min_mapq", action="store", default=2, type=int,\
        help='Minimum mapq score of EITHER read in a pair to use that pair. (>=) ')
    parser.add_argument("--min_insert", action="store", default=50, type=int,\
        help='Minimum insert size of two reads - default is 50 bp. True insert size includes read lengths: if two reads are 50bp each and overlap completely, their insert will be 50. (>=)')
    parser.add_argument("--max_insert", action="store", default=1500, type=float, \
        help='Maximum insert size between two reads - default is 1500 bp (<=)')
    parser.add_argument("--min_breadth_cov", action="store", default="0,0", type=str,\
        help='OPTIONAL: comma separated min breadth and a specific coverage required to generate output. E.g., 50,5 will require 50 percent of the genome to be at least 5x coverage to generate tables')

    # Parse
    if (len(args) == 0 or args[0] == '-h' or args[0] == '--help'):
        parser.print_help()
        sys.exit(0)
    else:
        return parser.parse_args(args)

if __name__ == '__main__':
    Controller().main(sys.argv[1:])
