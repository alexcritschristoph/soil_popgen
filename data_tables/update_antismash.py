scafs = set()
from Bio import SeqIO
import glob

for fn in glob.glob('../representative_genomes/*.fasta'):
	for scaf in SeqIO.parse(fn, "fasta"):
		scafs.add(scaf.id)

f = open('antismash_genes.tsv')
antismash = {}
print(f.readline())
for line in f.readlines():
	scaf = line.split()[3]
	if scaf in scafs:
		print(line.strip())
