# Bait sequences

Bromeliad1776 is a target sequencing probe set to enrich 1776 single and multi-copy genic regions across Bromeliaceae.


## The files
- align_and_trim.sh - Generic pipeline to align and trim sequenced captures DNA against the probe target file. Demultiplexing with deML, trimming with trim_galore and fastqc reports, aligned with bowtie2 - VCF called with freebayes.
- calculate_bait_target_specifity.sh - Pipeline to calculate target specifity, given a list of aligned and quality trimmed bams and a bed file with target coordinates, use samtools and bedtools intersect to look at % of reads mapped to target.
- pineapple.20150427.annotat.gff3 - gff of _Ananas_ genome used to design the set. Note: _Ananas comosus_ v.3 (available of phytozome)
- Table_S1_All_Genes_Details.csv - Details about the genes included in the final design (annotations follow those of *A. comosus*)
- legend_Table_S1.csv - legend for table S1
