# Soil bacterial populations are shaped by recombination and gene-specific selection across a grassland meadow
Reproducible scripts and notebooks for 2019 paper on population genetics in metagenomes

## Structure

*inStrain_lite* - A `pip` installable python program that will call SNPs and calculate nucleotide diversity, linkage, coverage+breadth, and Fst.
Meant to be reusable for other projects, eventually will be published as its own independent repository. 

*Data_tables* - Most of the important data tables used in this study.

*Notebooks* - notebooks used to generate each figure in this paper.

*Representative_genomes* - representative genome sequences for each species used in this study, along with predicted genes.

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

A BAM file for each species containing all filtered reads assigned to that species from the meadow was created using `combine_filter_bams.py`.  These BAMs were then used to analyze total nucleotide diversity and SNP linkage across the meadow. These BAMs were then profiled with the `inStrain_lite` script with the following parameters:
```
inStrain_lite -p 48 -s 30 -c 0.96 ../meadow_wide/14_0903_02_20cm_Proteobacteria_56_68_14_filtered_sort.bam  ../representative_genomes
/14_0903_02_20cm_Proteobacteria_56_68_14.fasta
```

### Per-sample population profiling

Each individual sample was profiled separately for the purpose of calculating nucleotide diversity per sample / plot, and Fst between plots. The inStrain commands run take the form: 
```
inStrain_lite -p 48 -s 30 -c 0.96 --min_breadth_cov 0.5,5 ../bams/all_14_0903_02_20cm_sorted.bam ../representative_genomes/14_0903_02_20cm_Proteobacteria_56_68_14.fasta
```
586 out of 1140 total (19 genomes*60 samples) genome+sample pairs passed the minimum requirement of at least 5x filtered read coverage across at least 50% of the representative genome to be included in these downstream analyses. 

Data from samples were then aggregated into various groupings (by replicate, plot, block, and all) using `make_merge_commands.py` which generates commands for all of the possible combinations to run `./inStrain_lite/combine_samples.py`. Merged samples were used for Figure 3 and for FST data (merged by block).

### Gene statistics

Genes were called on all representative genomes with `prodigal` - output is in `./representative_genomes/fna/*`. The script `gene_statistics.py` was run on each of the per-sample population profiles, each aggregated profile, and the meadow-wide profiles (commands found in scripts called `run_gene_commands.sh` in each respective `./data/` directory). For each gene, multi-allelic nucleotide diversity was reported, and point mutations were labeled as N or S using a naive point mutation test for a change in translation using the default genetic code.

### Organism abundances

Organism relative abundances were taken by counting the total number of reads that mapped to each genome in each sample (`./data_tables/bin_reads.txt`) and dividing by the total number of reads per sample (`bin_total_reads.tsv`). 

### Fst and Tajima's D

Fst between and Tajima's D within the three meadow blocks were calculated using the script `fst.py` for every gene. This script requires the package scikit-allel which it uses to calculate these statistics. Statistics were calculated on SNPs segregating within or between either population (sps were jointly called across both populations). Fst was calculated using the Hudson method as recommended by Bhatia et al. (2013).
