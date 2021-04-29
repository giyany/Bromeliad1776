#!/bin/bash
# given a list of bam filenames (aligned and quality trimmed) [-a] and a bed with target coordinates [-b], 
#use samtools and bedtools intersect to look at % of reads mapped to target

usage() { echo "Usage: $0 [-a <list of aligned bams files>] [-b <bed file of targets>]" 1>&2; exit 1; }


while getopts ":a:b:" opt; do
  case $opt in
    a)
      a=${OPTARG}
      ;;
     b)
	  b=${OPTARG}
      ;;
    \?)
      usage >&2
      ;;
  esac
done

if [ -z "${a}" ] || [ -z "${b}" ]; then
    usage
fi

 while read file; do
 	echo "generating stats on $file";
 	printf "$file" >> target_specifity_report2.txt; printf "\t" >> target_specifity_report2.txt;
 	samtools stats "$file" | awk 'NR==14 {printf $4}' >> target_specifity_report2.txt; printf "\t" >> target_specifity_report2.txt;
 	echo "intersecting $file with targets"
 	bedtools intersect -a "$file" -b "$b" > "$file".ontarget.bam;
 	echo "calculating specificity for $file";
 	samtools stats "$file".ontarget.bam | awk 'NR==14 {printf $4}' >> target_specifity_report2.txt; printf "\t" >> target_specifity_report2.txt;
 	printf '\n' >> target_specifity_report2.txt;
 	rm "$file".ontarget.bam;
 done < "$a"
 	awk '{print $1"\t"$2"\t"$3"\t"$3/$2}' target_specifity_report2.txt > target_specifity_report.txt;
 	rm target_specifity_report2.txt
