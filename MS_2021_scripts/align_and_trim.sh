#!/bin/bash

indlist=$1

#preparing file for report
printf "%s\t" "ind" "reads_paired" "aligned_reads" "mapped_paired" "AVG_qual" "uniquely_mapped" "surviving_QC" "duplication_per"> alignment_report.txt
echo " " >> alignment_report.txt;

#trim
mkdir reports;
while read ind; do
         fastqc -t 12 -f fastq "$ind"_R1.fq "$ind"_R2.fq;
         java -jar ~/Programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 12 -phred33 -trimlog "$ind"_trim_log "$ind"_R1.fq "$ind"_R2.fq "$ind"_R1_paired.fq "$ind"_R1_unpaired.fq "$ind"_R2_paired.fq "$ind"_R2_unpaired.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50;
         mv "$ind"_trim_log reports;
         fastqc -t 12 "$ind"_R1_paired.fq "$ind"_R2_paired.fq;
 done < $indlist

#awkward code for collecting info about read trimming
 while read ind; do
         unzip "$ind"_R1_fastqc.zip;
         echo "R1 before trim" >> reports/"$ind"_Trim_Sequences_Count.txt;
         less "$ind"_R1_fastqc/fastqc_data.txt | grep "Total Sequences" >> reports/"$ind"_Trim_Sequences_Count.txt;
         rm -r "$ind"_R1_fastqc;
         mv "$ind"_R1_fastqc.zip reports;
         unzip "$ind"_R2_fastqc.zip;
         echo "R2 before trim" >> reports/"$ind"_Trim_Sequences_Count.txt;
         less "$ind"_R2_fastqc/fastqc_data.txt | grep "Total Sequences" >> reports/"$ind"_Trim_Sequences_Count.txt;
         rm -r "$ind"_R2_fastqc;
         mv "$ind"_R2_fastqc.zip reports;
         unzip "$ind"_R1_paired_fastqc.zip;
         echo "R1 after trim" >> reports/"$ind"_Trim_Sequences_Count.txt;
         less "$ind"_R1_paired_fastqc/fastqc_data.txt | grep "Total Sequences" >> reports/"$ind"_Trim_Sequences_Count.txt;
         rm -r "$ind"_R1_paired_fastqc;
         mv "$ind"_R1_paired_fastqc.zip reports;
         unzip "$ind"_R2_paired_fastqc.zip;
         echo "R2 after trim" >> reports/"$ind"_Trim_Sequences_Count.txt;
         less "$ind"_R2_paired_fastqc/fastqc_data.txt | grep "Total Sequences" >> reports/"$ind"_Trim_Sequences_Count.txt;      
         rm -r "$ind"_R2_paired_fastqc;
         mv "$ind"_R2_paired_fastqc.zip reports;
         mv "$ind"_R1_fastqc.html reports;
         mv "$ind"_R2_fastqc.html reports;
         mv "$ind"_R1_paired_fastqc.html reports;
         mv "$ind"_R2_paired_fastqc.html reports;
 done < $indlist

#aligning to reference
while read ind; do
	do echo "bowtie aligning $ind";
	bowtie2 --very-sensitive-local -x /media/computer/097c7584-78ef-4295-907e-df3c187190c3/Science/Data/Refgenome/Tfas_bowtieindex -1 "$ind"_R1_paired.fq -2 "$ind"_R2_paired.fq -S "$ind"_aligned_Tfasc.sam -p 14;
	#collect stats on alignment
	samtools stats "$ind"_aligned_Tfasc.sam > "$ind"_aligned_Tfasc_sam_stats.txt;
	#convert to bam, keep uniquely mapped reads only, filter by mq>10 & sort
	samtools view -h -b -q 10 "$ind"_aligned_Tfasc.sam | samtools sort -o "$ind"_asmq10.bam;
	#add read group info. This is specific to my data, other users would want to modify
	java -jar ~/Programs/picard.jar AddOrReplaceReadGroups I="$ind"_asmq10.bam o="$ind"_asmq10rg.bam RGLB=WGD RGPL=illumina RGPU=Lib1 RGSM="$ind" RGID="$ind";
	#mark duplicates. notice duplicates here are MARKED NOT REMOVED
	java -jar ~/Programs/picard.jar MarkDuplicates I="$ind"_asmq10rg.bam o="$ind"_asmq10rgd.bam M="$ind"_dup_metrics.txt;
	#redirect info about alignment into txt report
	printf "%s\t" "$ind" >> alignment_report.txt;
	awk 'BEGIN {ORS="\t"} NR==10 {print $3} NR==14 {print $4} NR==15 {print $6} NR==32 {print $4}' reports/"$ind"_aligned_Tfasc_sam_stats.txt >> alignment_report.txt;
	samtools view -F 4 "$ind"_aligned_Tfasc.sam | grep -v "XS:" | wc -l >> alignment_report.txt;
	truncate -s -1 some_text.txt; printf "%s\t" "" >> alignment_report.txt;
	awk 'BEGIN {ORS="\t"} NR==8 {print $3*2} NR==8 {print $9}' "$ind"_dup_metrics.txt >> alignment_report.txt;
	echo " " >> alignment_report.txt;
	mv "$ind"_aligned_Tfasc_sam_stats.txt reports;
	mv "$ind"_dup_metrics.txt reports;
done < $indlist

# Call variants with freebayes, example (code for cluster computing)

TEMPDIR=/gpfs/data/fs71400/yardeni/Angio353/vcfcall_Ananas/tempdir freebayes-parallel < /gpfs/data/Ananas_resources/Ananas.target.regions.10k.by.cov.LG.only 32 \ --report-monomorphic --use-best-n-alleles 120 --limit-coverage 600 -g 50000 -f pineapple.20150427.fasta --populations popfile.txt *.bam > Angio353vsAnanas_wmonosites_raw.vcf
