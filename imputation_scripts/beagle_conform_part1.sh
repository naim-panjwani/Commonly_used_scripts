#!/bin/bash
# Syntax: bash beagle_conform_part1.sh <vcf_sample_filename w/out .vcf.gz> <output_filename_prefix>
# Example: Suppose we have sample files 07_UKRE.chr1.vcf.gz, *.chr2.*, etc
#          bash beagle_conform_part1.sh 07_UKRE 08_UKRE

ref_dir="1000Genomes_beagle_v4_format/"
sample=$1

for chr in $(seq 1 11); do
  echo "chr${chr}"
  ref="${ref_dir}chr${chr}.1kg.ref.phase1_release_v3.20101123.vcf.gz"
  samplename="${sample}.chr${chr}.vcf.gz"
  conformed="$2.chr${chr}.conformed"
# out="${sample}.chr${chr}.imputed"
  
  nohup java -jar conform-gt.r1174.jar ref=${ref} gt=${samplename} chrom=${chr} out=${conformed} >conform_chr${chr}.out &
done

