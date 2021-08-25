#!/bin/bash

indlist=$1

deML -i P5P7_for_demulti --bamtags BC,QT,B2,Q2 -o Angio353-data-demulti_deML.bam -s demult_stats.txt -e demult_unassigned.txt Angio353-data_raw.bam

#preparing file for report
printf "%s\t" "ind" "reads_paired" "aligned_reads" "mapped_paired" "AVG_qual" "uniquely_mapped" "surviving_QC" "duplication_per"> alignment_report.txt
echo " " >> alignment_report.txt;

while read ind; do
	bamtools filter -tag RG:"$ind" -in Angio353-data-demulti_deML.bam  -out "$ind"_Ang353.bam;
	samtools stats "$ind"_Ang353.bam > "$ind"_raw_stats.txt;
	bedtools bamtofastq -i "$ind"_Ang353.bam -fq "$ind"_Ang353_R1.fq -fq2 "$ind"_Ang353_R2.fq;
done < ind_list.txt

#trim and report
mkdir reports;
while read ind; do
	fastqc -t 12 -f fastq "$ind"_R1.fq "$ind"_R2.fq;
	trim_galore --fastqc --cores 8 --retain_unpaired --paired "$ind"_R1.fq "$ind"_R2.fq;
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
	bowtie2 --very-sensitive-local -x /media/computer/097c7584-78ef-4295-907e-df3c187190c3/Science/Data/Refgenome/Acom_bowtieindex -1 "$ind"_R1_paired.fq -2 "$ind"_R2_paired.fq -S "$ind"_aligned_Acom.sam -p 14;
	#collect stats on alignment
	samtools stats "$ind"_aligned_Acom.sam > "$ind"_aligned_Acom_sam_stats.txt;
	#convert to bam, keep uniquely mapped reads only, filter by mq>10 & sort
	samtools view -h -b -q 10 "$ind"_aligned_Acom.sam | samtools sort -o "$ind"_asmq10.bam;
	#add read group info. This is specific to my data, other users would want to modify
	java -jar ~/Programs/picard.jar AddOrReplaceReadGroups I="$ind"_asmq10.bam o="$ind"_asmq10rg.bam RGLB=WGD RGPL=illumina RGPU=Angio353 RGSM="$ind" RGID="$ind";
	#mark duplicates. notice duplicates here are MARKED NOT REMOVED
	java -jar ~/Programs/picard.jar MarkDuplicates I="$ind"_asmq10rg.bam o="$ind"_asmq10rgd.bam M="$ind"_dup_metrics.txt;
	#redirect info about alignment into txt report
	printf "%s\t" "$ind" >> alignment_report.txt;
	awk 'BEGIN {ORS="\t"} NR==10 {print $3} NR==14 {print $4} NR==15 {print $6} NR==32 {print $4}' reports/"$ind"_aligned_Acom_sam_stats.txt >> alignment_report.txt;
	samtools view -F 4 "$ind"_aligned_Acom.sam | grep -v "XS:" | wc -l >> alignment_report.txt;
	truncate -s -1 some_text.txt; printf "%s\t" "" >> alignment_report.txt;
	awk 'BEGIN {ORS="\t"} NR==8 {print $3*2} NR==8 {print $9}' "$ind"_dup_metrics.txt >> alignment_report.txt;
	echo " " >> alignment_report.txt;
	mv "$ind"_aligned_Acom_sam_stats.txt reports;
	mv "$ind"_dup_metrics.txt reports;
done < $indlist

# Call variants with freebayes, example (code for cluster computing)

TEMPDIR=/gpfs/data/fs71400/yardeni/Angio353/vcfcall_Ananas/tempdir freebayes-parallel < /gpfs/data/Ananas_resources/Ananas.target.regions.10k.by.cov.LG.only 32 \ --report-monomorphic --use-best-n-alleles 120 --limit-coverage 600 -g 50000 -f pineapple.20150427.fasta --populations popfile.txt *.bam > Angio353vsAnanas_wmonosites_raw.vcf

#some filtering steps
#break MNPs down and fill in AC, AN tag                                                                                                                                                                      
                                                                                                                                                                                                             
vcfallelicprimitives -kg Angio353_Ananas_aligned_allsamples_vcf_raw_ALL_LGs.vcf.bgzip > Angio353_Ananas_aligned_allsamples_vcf_raw_primitives.vcf;                                     
bcftools plugin fill-AN-AC --threads 16 Angio353_Ananas_aligned_allsamples_vcf_raw_primitives.vcf | bcftools sort - | bgzip > Angio353_Ananas_aligned_allsamples_vcf_raw_primitives_ANAC.vcf.gz              
                                                                                                                                                                                                             
#remove indels and SNPs 3bp around indels, allow max missing 40% in site                                                                                                                                     
                                                                                                                                                                                                             
bcftools filter --SnpGap 3 -i 'F_MISSING<0.4' Angio353_Ananas_aligned_allsamples_vcf_raw_primitives_ANAC.vcf.gz | bcftools view --exclude-types indels | bgzip > Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4.vcf.gz; tabix Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4.vcf.gz;                                                                                      
                                                                                                                                                                                                             
#depth for genotype > 10 for site > 100                                                                                                                                                                      
                                                                                                                                                                                                             
~/Programs/vcflib/bin/vcffilter -f "DP > 500" -g "DP > 10" Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4.vcf.gz > Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500.vcf;                                                                                                                                                                                           
                                                                                                                                                                                                             
#at least 10 copies of the alternative allele present [but keep monomorphic sites]. This step removes all non-variant sites, so I extract them separatedly then concatenated together                         
                                                                                                                                                                                                             
bcftools filter -i 'AC>9 | AC==0' Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500.vcf | bcftools sort - | bgzip > Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_Nalt.vcf.gz; tabix Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_Nalt.vcf.gz;                                                                         
bcftools view -C0 Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500.vcf | bcftools sort - | bgzip > Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_monoonly.vcf.gz; tabix Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_monoonly.vcf.gz;                                                                                 
                                                                                                                                                                                                             
bcftools concat -a -D Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_Nalt.vcf.gz Angio353_Ananas_aligned_allsamples_vcf_filter_SnpGap_indels_missing0.4_DP10_500_monoonly.vcf.gz | bcftools sort - | bgzip > Angio353_Ananas_aligned_allsamples_withmono_filter_SnpGap_indels_missing0.4_DP10_500_Nalt.vcf;
