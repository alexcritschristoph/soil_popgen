import numpy as np 
import pandas as pd 
from tqdm import tqdm
from inStrain_lite import call_snv_site
from inStrain_lite import SNPprofile

def calc_fst_weir()
def calc_fst(allele_freq1, allele_freq2):
    Hs = allele_freq1*(1-allele_freq1) + allele_freq2 * (1-allele_freq2)
    avg_freq = np.mean([allele_freq1, allele_freq2])
    Ht = 2*avg_freq*(1-avg_freq)
    Fst = (Ht - Hs) / Ht
    return(Fst)
# def FST_H_pairwise (col):
#     NpopA = float(col[0])
#     NpopB = float(col[2])
#     popAcount= int(col[1])
#     popBcount= int(col[3])  
#     pop1freq = popAcount / float(NpopA )
#     pop2freq = popBcount / float(NpopB )
#     Npop1 = NpopA
#     Npop2 = NpopB
#     T_1=(pop1freq-pop2freq)**2 - ((pop1freq*(1.0-pop1freq))/(Npop1-1)) - ((pop2freq*(1.0-pop2freq))/(Npop2-1))
#     T_2=(pop1freq*(1.0-pop2freq)) + (pop1freq*(1.0-pop2freq))
#     return (T_1,T_2)


# [docs]
# def hudson_fst(ac1, ac2, fill=np.nan):

#     # calculate these once only
#     an1 = np.sum(ac1, axis=1)
#     an2 = np.sum(ac2, axis=1)

#     # calculate average diversity (a.k.a. heterozygosity) within each
#     # population
#     within = (mean_pairwise_difference(ac1, an1, fill=fill) +
#               mean_pairwise_difference(ac2, an2, fill=fill)) / 2

#     # calculate divergence (a.k.a. heterozygosity) between each population
#     between = mean_pairwise_difference_between(ac1, ac2, an1, an2, fill=fill)

#     # define numerator and denominator for Fst calculations
#     num = between - within
#     den = between

#     return num, den

def main(args):

    ## Step 0: get null model for SNP calling
    null_loc = os.path.dirname(__file__) + '/helper_files/combined_null1000000.txt'
    null_model = generate_snp_model(null_loc)
    P2C = {'A':0, 'C':1, 'T':2, 'G':3}
    C2P = {0:'A', 1:'C', 2:'T', 3:'G'}


    ## Step 1: build new counts table from all objects
    s_final = SNVprofile()
    s_final.filename = args.output
    i = 0
    counts_per_block = {}
    s1 = SNVprofile()
    print("loading " + args.input[0])
    s1.load(args.input[0])

    s_final.scaffold_list = s1.scaffold_list
    s_final.counts_table = s1.counts_table
    
    s2 = SNVprofile()
    print("loading " + args.input[1])
    s2.ln(args.input[1])

    for scaf in s2.scaffold_list:
        if scaf not in s_final.scaffold_list:
            sys.exit("Error: scaffold " + scaf + " in " + fn + " not found in initial file. Your inStrain objects were probably not run on the same FASTA.")

    scaf_counter = 0
    for scaf in s1.counts_table:
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
                    snps_table['scaffold'] = s_final.scaffold_list[scaf_counter]
                    snps_table['position'] = pos_counter
                    snps_table['varBase'] = snp
                    snps_table['conBase'] = 
                    snps_table['var_freq_pop1']
                    snps_table['var_freq_pop2']

            pos_counter += 1 # 0 based positions!!
        scaf_counter += 1

    # Step 3: Save new FST_SNP table to disk.
    SNPTable = pd.DataFrame(snp_table)
    SNPTable.to_csv(args.output + '.SNVs.tsv', index=False, sep='\t')
    
    # Step 4: Collate FST by gene

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description= """Calculates Fst """, formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument('input', nargs='+', help="two inStrain object prefixes, space separated. (don't include `.`, only prefixes - must be in the local directory).")
    parser.add_argument("-o", "--output", action="store", required=True, \
        help='Output prefix ')

    args = parser.parse_args()
    main(args)