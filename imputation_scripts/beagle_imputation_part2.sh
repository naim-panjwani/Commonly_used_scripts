#!/bin/sh
# Syntax: bash beagle_imputation_part2.sh <sample_filename w/out .vcf.gz> <output filename prefix> <starting chromosome num; if chrX, specify 23> <num simultaneous processes to run>
# Example: Suppose we have UKRE.chr1.conformed.vcf.gz, *.chr2.*, etc. files, the last chromosome that was run by part1 was chr13, and most of the resources on Chuck
#          nohup bash beagle_imputation_part2.sh 09_UKRE 10_UKRE 14 13 >imputation_part2.out &
# NOTE: YOU MAY START PART2 AT ANYTIME; IT WILL NOT EXCEED THE NUMBER OF JAVA PROCESSES SPECIFIED

count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`
sample=$1
#echo "Starting chromosome to impute:"
chr=$3
samplename="${sample}.vcf.gz"
ref_dir="1000Genomes_beagle_v4_format/"
max_num_processes=$4

while [[ $chr -le 23 ]]; do
  while [[ $count -ge ${max_num_processes} ]]; do 
    sleep 60
    count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`
  done

  ref="${ref_dir}chr${chr}.1kg.ref.phase1_release_v3.20101123.vcf.gz"
  conformed="${sample}.chr${chr}.conformed"
  out="$2.chr${chr}.imputed"  

  if [[ $chr -eq 23 ]] ; then
#    nohup java -Xmx30000m -jar beagle.r1398.jar ref=${ref_dir}chrX.1kg.ref.phase1_release_v3.20101123.vcf.gz gt=${sample}.chrX.conformed.vcf.gz out=${sample}.chrX.imputed chrom=X >impute_chrX.out &
    nohup java -jar beagle.r1398.jar ref=${ref_dir}chrX.1kg.ref.phase1_release_v3.20101123.vcf.gz gt=${sample}.chrX.conformed.vcf.gz out=$2.chrX.imputed chrom=X >impute_chrX.out &
    echo "Chr X started"
    sleep 5
  else
#    nohup java -Xmx30000m -jar beagle.r1398.jar ref=${ref} gt=${conformed}.vcf.gz out=${out} chrom=${chr} >impute_chr${chr}.out &
    nohup java -jar beagle.r1398.jar ref=${ref} gt=${conformed}.vcf.gz out=${out} chrom=${chr} >impute_chr${chr}.out &
    echo "Chr ${chr} started"
    sleep 5
  fi
  chr=$[$chr + 1]
  count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`
done
echo Done sending imputation for last chromosome >DONE_job_sending
