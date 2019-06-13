import os
import sys
import pysam
import logging
import argparse
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from inStrain_lite.utilities import setup_logger
from collections import defaultdict

def get_fasta(fasta_file = None):
    scaffolds = []
    total_length = 0
    for rec in SeqIO.parse(fasta_file, "fasta"):
        start = 0
        total_length += len(rec.seq)
        scaffolds.append(str(rec.id))
    return scaffolds, total_length

def filter_paired_reads(bams, fasta, min_ani = 0.97, min_mapq = 2, min_insert = 50, max_insert = 1500, write_bam = False):
    '''
    Filter reads from a .bam file
    Returns:
        pair2info - dictionary of read pair -> (mismatches, insert distance, mapq score, combined length)
    '''

    logging.info("Copying header for new bam...")
    samfile0 = pysam.AlignmentFile(bams[0])
    samfile_out = pysam.AlignmentFile(fasta.split("/")[-1].split(".")[0] + "_filtered.bam", "wb", template=samfile0)
    samfile0.close()
    logging.info("getting fasta")
    scaffolds, fasta_length = get_fasta(fasta)
    total = mapped_pairs = mapq_good = insert_good = 0
    insert_sizes = []
    read_lengths = []
    filtered = set() # Information on pairs

    for bam in bams:
        samfile = pysam.AlignmentFile(bam)

        for scaff in tqdm(scaffolds, desc='Filtering Reads'):
            read_data = {} # Information on the first pair of each read
            for read in samfile.fetch(scaff):
                total += 1
                # If we've seen this read's pair before
                if read.query_name in read_data:
                    # Make sure that the pair is on the same scaffold and that it's mapped at all
                    if ((read_data[read.query_name]['scaf'] == scaff) & (read.get_reference_positions() != [])):
                        mapped_pairs += 1
                        pairMM = float(read_data[read.query_name]['read'].get_tag('NM')) + float(read.get_tag('NM')) #number of mismatches in pair
                        mapped_read_lengths = float(read_data[read.query_name]['read'].infer_query_length() + read.infer_query_length()) #total length of pair
                        if read.get_reference_positions()[-1] > read_data[read.query_name]['read'].get_reference_positions()[0]:
                            pair_inserts = read.get_reference_positions()[-1] - read_data[read.query_name]['read'].get_reference_positions()[0] #insert distance
                        else:
                            pair_inserts = read_data[read.query_name]['read'].get_reference_positions()[-1] - read.get_reference_positions()[0] #insert distance
                        pair_mapq = max(read.mapping_quality, read_data[read.query_name]['read'].mapping_quality) #pair mapq
                        pair_ani =  1 - (pairMM / mapped_read_lengths) #pair %ANI to reference
                        insert_sizes.append(pair_inserts)
                        read_lengths.append(mapped_read_lengths)
                        # Final filter
                        if pair_inserts >= min_insert and pair_inserts <= max_insert:
                            insert_good += 1
                        if pair_mapq >= min_mapq:
                            mapq_good += 1
                            if pair_ani >= min_ani:
                                filtered.add(read.query_name)
                                samfile_out.write(read_data[read.query_name]['read'])
                                samfile_out.write(read)

                # Add this read, in future search for its mate
                elif read.get_reference_positions() != []: # don't use unmapped reads:
                    read_data[read.query_name] = {"read": read, "scaf": scaff}
        total = total / 2
        logging.info("Total expected number of read pairs\t" + str(total))
        logging.info("Mean read pair sequences length\t" + str(np.mean(read_lengths)))
        logging.info("Total FASTA length\t" + str(fasta_length))
        logging.info("Expected total coverage\t" + str(float(total)*np.mean(read_lengths) / fasta_length))
        logging.info("Mapped read pairs\t" + str(mapped_pairs) + "\t" + str(int(100*mapped_pairs / total)) + "%")
        logging.info("Median end-to-end insert length\t" + str(np.median(insert_sizes)))
        logging.info("Read pairs which pass insert distance filters:\t" + str(insert_good) + "\t" + str(int(100*float(insert_good) / total)) + "%")
        logging.info("Read pairs which also meet min_mapq of " + str(min_mapq) + "\t" + str(mapq_good) + "\t" + str(int(100*float(mapq_good) / total)) +  "%")
        logging.info("Read pairs which also pass final read pair PID >" + str(min_ani) + "%\t" + str(len(filtered)) + "\t" + str(int(100*len(filtered) / total)) + "%")
        logging.info("Final expected coverage\t" + str(float(len(filtered)) * np.mean(read_lengths) / fasta_length))
        samfile.close()
    
    samfile_out.close()
    logging.info("sorting new bam")
    pysam.sort("-o", fasta.split("/")[-1].split(".")[0] + "_filtered_sort.bam", fasta.split("/")[-1].split(".")[0] + "_filtered.bam")
    os.system('rm ' + fasta.split("/")[-1].split(".")[0] + "_filtered.bam")
    return filtered

if __name__ == '__main__':
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser(description="""Filters multiple BAM files and generates a combined BAM of the filtered reads from each BAM.""", formatter_class=argparse.RawTextHelpFormatter)

    # Required positional arguments
    parser.add_argument("bams", nargs='+', help="list of sorted bams files")
    parser.add_argument("fasta", help="Fasta file the bam is mapped to")
    # Optional arguments
    parser.add_argument("-m", "--mismatch_threshold", action="store", default=0.97, \
        help='Minimum percent identity of read pairs to consensus to use the reads - default is to run at 0.97.')
    parser.add_argument("-q", "--min_mapq", action="store", default=2, \
        help='Minimum required mapq score of EITHER read in a pair to use that pair. Default: 2.')
    parser.add_argument("-l", "--max_insert", action="store", default=1500, \
        help='Maximum insert size between two reads - default is 1500 bp.')
    parser.add_argument("-u", "--min_insert_length", action="store", default=50, \
        help='Minimum insert size between two reads - default is 50 bp.')
    parser.add_argument("-g", "--generate_bam", action="store_true", default=True, \
        help='Include to create a new filtered BAM.')
    parser.add_argument('--log', action='store', default='filter.log', \
        help ="File to log results to.")

    args = parser.parse_args()
    setup_logger(args.log)
    filter_paired_reads(args.bams, args.fasta, float(args.mismatch_threshold), int(args.min_mapq), int(args.min_insert_length), int(args.max_insert), write_bam=args.generate_bam)

