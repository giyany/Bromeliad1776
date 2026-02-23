# Scripts for analysis in 2021 Manuscript

These scripts assist to replicate the analysis described in Yardeni et al. 2021 (https://doi.org/10.1111/1755-0998.13523).

## The files
- align_and_trim.sh - Generic pipeline to align and trim sequenced captures DNA against the probe target file. Demultiplexing with deML, trimming with trim_galore and fastqc reports, aligned with bowtie2 - VCF called with freebayes.
- calculate_bait_target_specifity.sh - Pipeline to calculate target specifity, given a list of aligned and quality trimmed bams and a bed file with target coordinates, use samtools and bedtools intersect to look at % of reads mapped to target.
- vcf2genocountsmatrix.py - A script for generating summary statistics of the final SNP sets, specifically the total number of SNPs, the proportion of on-target SNPs and the proportion of SNPs in specific genomic contexts 
