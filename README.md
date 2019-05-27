# soil_popgen
Reproducible scripts and notebooks for 2019 paper on population genetics in metagenomes

## Structure

*inStrain_lite* - A `pip` installable python program that will call SNPs and calculate nucleotide diversity, linkage, coverage+breadth, and Fst.
Meant to be reusable for other projects, eventually will be published as its own independent repository. 

*Data* - contains some of the usable data files for this paper.

*Notebooks* - notebooks used to generate each figure in this paper.

## Reproducibility

### Genome dereplication
dRep was run on all 10,538 genome bins obtained, with the secondary clustering threshold `-sa 0.97`.
CheckM genome information for all bins is shown in `genomeInfo.csv`.

### Pan-genome decontamination 
The representative genomes used in this study were decontaminated using the assembled replicate genomes for each species population. Roary was run with default settings on the set of genomes for each species. Contigs with at least 50% of their genes found in no more than 25% of each genome set were discarded as potential contamination. 

### Read filtering

Read filtering is performed using the custom script `filter_reads.py`. Briefly, a bowtie2 index is created using representative members for all 664 dereplicated genomes. Reads were then mapped to this index using bowtie2 default settings with the expected insert distance `-X 1000`. The resulting BAM files are then filtered for reads that meet the following criteria:
(1) Both read pairs map to the same scaffold within 1500 bp of each other (end to end insert size), with a percent identity of 96%.
(2) At least one of the read pairs has a mapq score > 1, e.g. this is a uniquely best mapping for this read (and therefore probably for this pair) in the index.

Reads in this dataset are on average ~200 bp. (HiSeq 2500 sequencing that aimed for 2x250 bp reads, but the trailing ends of the reads were often low quality and so they were usually trimmed to ~200 bp. 

### Meadow-wide population profiling.

A BAM file for each species containing all filtered reads assigned to that species from the meadow was created using `combine_filter_bams.py` in `./meadow_wide/`.  These BAMs were then used to analyze SNP linkage across the meadow. 

### Per-sample population profiling


