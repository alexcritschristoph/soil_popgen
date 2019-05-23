import argparse
from inStrain import SNPprofile
from inStrain import call_snv_site
from inStrain import generate_snp_model
import sys
import os
import pandas as pd
from collections import defaultdict
from tqdm import tqdm
import numpy as np

def main(args):

    ## Step 0: get null model for SNP calling
    null_loc = os.path.dirname(__file__) + './inStrain/helper_files/combined_null1000000.txt'
    null_model = generate_snp_model(null_loc)
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}


    ## Step 1: build new counts table from all objects
    s_final = SNVprofile()
    s_final.filename = args.output
    i = 0 
    for fn in args.input:
        s = SNVprofile()
        print("loading " + fn)
        s.load(fn)

        if i == 0: # first iteration through loop, first file provides the base objects
            s_final.scaffold_list = s.scaffold_list
            s_final.counts_table = s.counts_table
        else:
            for scaf in s.scaffold_list:
                if scaf not in s_final.scaffold_list:
                    sys.exit("Error: scaffold " + scaf + " in " + fn + " not found in initial file. Your inStrain objects were probably not run on the same FASTA.")
            
            scaf_counter = 0
            for scaf in s.counts_table:
                s_final.counts_table[scaf_counter] += scaf
                scaf_counter += 1
        i += 1


    # Step 2: call all SNPs for new object
    snp_table = defaultdict(list)
    scaf_counter = 0
    pos_counter = 0
    for scaf in tqdm(s_final.counts_table, desc='Calling new SNVs...'):
        for counts in scaf:
            snp = call_snv_site(counts, min_cov = args.min_coverage, min_freq = args.min_freq, model=null_model)

            if snp: # means that there was coverage at this position
                if snp != -1: # means this is a SNP
                    # calculate varBase
                    scaf = s_final.scaffold_list[scaf_counter]
                    
                    snp, varbase = major_minor_allele(counts)
                    
                    snp_table['scaffold'].append(scaf)
                    snp_table['position'].append(pos_counter)
                    for b, c in zip(['A', 'C', 'T', 'G'], counts):
                        snp_table[b].append(c)
                    snp_table['conBase'].append(snp)
                    snp_table['varBase'].append(varbase)
                    baseCoverage = np.sum(counts)
                    snp_table['baseCoverage'].append(baseCoverage)
                    snp_table['conFreq'].append(float(counts[P2C[snp]]) / baseCoverage)
                    snp_table['varFreq'].append(float(counts[P2C[varbase]]) / baseCoverage)
                   
            pos_counter += 1 # 0 based positions!!
        scaf_counter += 1

    # Step 3: Save new object and new SNP table to disk.
    SNPTable = pd.DataFrame(snp_table)
    SNPTable.to_csv(args.output + '.SNVs.tsv', index=False, sep='\t')
    s_final.store()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """
        Combines multiple inStrain objects for joint-sample snp calling and nucleotide diversity calculations.\n
        Cannot do combined linkage calculations, only joint-sample snp and nucleotide calculations.\n
        Returns a counts-only SNVprofile object along with a SNP table. 
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('input', nargs='+', help="multiple inStrain prefixes, space separated. (don't include `.`, only prefixes - must be in the local directory.")
    parser.add_argument("-o", "--output", action="store", required=True, \
        help='Output prefix ')

    parser.add_argument("-c", "--min_coverage", action="store", default=5, \
        help='Minimum SNV coverage')

    parser.add_argument("-f", "--min_freq", action="store", default=0.05, \
        help='Minimum SNP frequency to confirm a SNV (both this AND the 0.001 percent FDR snp count cutoff must be true).')

    args = parser.parse_args()

    main(args)

