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
