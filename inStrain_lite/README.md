# strains_analysis
Code for analyzing population genomics in genome-resolved metagenomes

Requires: pysam, tqdm, BioPython.

### Usage

```python strainRep2.py -s 5 -f 0.05 -l 0.96 sorted_indexed_bam.bam scaffolds.fasta -o  output --log  log.txt```

```python strainRep2.py -h``` 
actually is pretty helpful, that's all of the documentation.

output: 3 tables (and a big python object). Linkage table (showing snp linkage), frequency table (showing SNPs and their frequencies), and clonality table (showing the clonality and coverage of each position - from this gene clonality can be calculated and compared to the genome average) (edited) 
`-s 5` requires 5 reads to confirm a SNP, you can adjust depending on your coverage. `-f` means minimum snp frequency of 5%, `-l 0.96` means that read pairs must be 96% ID to reference. the statistics reported in the log file are also super useful
