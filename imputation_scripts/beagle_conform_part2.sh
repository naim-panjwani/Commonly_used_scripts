#!/bin/bash
# Syntax: bash beagle_conform_part2.sh <vcf_sample_filename w/out .vcf.gz> <output_filename_prefix> <starting chromosome> <num simultaneous processes to run>
# Example: Suppose we have sample files 07_UKRE.chr1.vcf.gz, *.chr2.*, etc
#          nohup bash beagle_conform_part2.sh 07_UKRE 08_UKRE 12 12 >beagle_conform_part2.out &
# NOTE: YOU MAY START PART2 AT ANYTIME; IT WILL NOT EXCEED THE NUMBER OF JAVA PROCESSES SPECIFIED

ref_dir="1000Genomes_beagle_v4_format/"
sample=$1
chr=$3
count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`

while [[ $chr -le 23 ]]; do
  while [[ $count -ge $4 ]]; do 
    sleep 60
    count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`
  done

  echo "chr${chr}"
  ref="${ref_dir}chr${chr}.1kg.ref.phase1_release_v3.20101123.vcf.gz"
  if [[ $chr -eq 23 ]] ; then
    samplename="${sample}.chrX.vcf.gz"
    conformed=$2.chrX.conformed
    nohup java -jar conform-gt.r1174.jar ref=${ref_dir}chrX.1kg.ref.phase1_release_v3.20101123.vcf.gz gt=${samplename} chrom=X out=${conformed} >conform_chrX.out &
  else
    samplename="${sample}.chr${chr}.vcf.gz"
    conformed="$2.chr${chr}.conformed"
  # out="${sample}.chr${chr}.imputed"
  
    nohup java -jar conform-gt.r1174.jar ref=${ref} gt=${samplename} chrom=${chr} out=${conformed} >conform_chr${chr}.out &
  fi
  chr=$[$chr + 1]
  count=`ps -C java |awk '{print $1}' |sed '1d' |wc -l`
done
echo "Done sending conform-gt jobs" >DONE_conform-gt_job_sending
