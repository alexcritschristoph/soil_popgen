### FUNCTIONS USED BY THE PROGRAM

from collections import defaultdict
import numpy as np
import logging

P2C = {'A':0, 'C':1, 'T':2, 'G':3}
C2P = {0:'A', 1:'C', 2:'T', 3:'G'}


def major_minor_allele(counts):
    d = {'A': counts[0], 'C': counts[1], 'T': counts[2], 'G': counts[3] }
    l = sorted(d, key=d.get, reverse=True)
    return l[0], l[1]

def generate_snp_model(model_file):
    f = open(model_file)
    model = defaultdict(int)
    for line in f.readlines():
        counts = line.split()[1:]
        coverage = line.split()[0]
        i = 0
        for count in counts:
            i += 1
            if float(count) < 1e-6:
                model[int(coverage)] = i
                break

    model.default_factory = lambda:max(model.values())
    return model

def get_base_counts(pileupcolumn, filtered_reads):
    '''
    From a pileupcolumn object, return a list with the counts of [A, C, T, G]
    '''
    counts = np.zeros(4, dtype=int)
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip and pileupread.alignment.query_name in filtered_reads:
            try:
                counts[P2C[pileupread.alignment.query_sequence[pileupread.query_position]]] += 1
            except KeyError: # This would be like an N or something not A/C/T/G
                pass
    return counts

def call_snv_site(counts, min_cov=5, min_freq=0.05, model=None):
    '''
    Determines whether a site has a variant based on its nucleotide count frequencies.
    Return:
        Base if SNP
        -1 if not SNP
        None if unCounted
    '''
    if model: #alexcc - so you can call this function from outside of the file
        model_to_use = model
    else:
        model_to_use = null_model

    total = sum(counts)
    if total >= min_cov:
        i = 0
        for c in counts:
            if c >= model_to_use[total] and float(c) / total >= min_freq:
                i += 1
        if i > 1:
            return C2P[np.argmax(counts)]
        else:
            return None
    else:
        return None

def setup_logger(loc):
    ''' set up logger such that DEBUG goes only to file, rest go to file and console '''

    # set up logging everything to file
    logging.basicConfig(level=logging.DEBUG,
                       format='%(asctime)s %(levelname)-8s %(message)s',
                       datefmt='%m-%d %H:%M',
                       filemode='w',
                       filename=loc)

    # set up logging of INFO or higher to sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    logging.debug("!"*80)
    logging.debug("Command to run was: {0}\n".format(' '.join(sys.argv)))
