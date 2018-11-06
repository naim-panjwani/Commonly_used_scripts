#!/bin/bash
# Syntax: bash beagle_imputation_part1.sh <sample_filename w/out .chr#.conformed.vcf.gz> <output_filename_prefix>
# Example: Suppose we have files that have been conformed-gt'd already: 09_UKRE.chr1.conformed.vcf.gz, *chr2*, etc.
#          bash beagle_imputation_part1.sh 09_UKRE 10_UKRE

sample=$1
samplename="${sample}.vcf.gz"
ref_dir="1000Genomes_beagle_v4_format/"

for chr in $(seq 1 11); do
  echo "chr${chr}"
  ref="${ref_dir}chr${chr}.1kg.ref.phase1_release_v3.20101123.vcf.gz"
  conformed="${sample}.chr${chr}.conformed"
  out="$22.chr${chr}.imputed"
  
#  nohup java -Xmx30000m -jar beagle.r1398.jar ref=${ref} gt=${conformed}.vcf.gz out=${out} chrom=${chr} >impute_chr${chr}.out &
  nohup java -jar beagle.r1398.jar ref=${ref} gt=${conformed}.vcf.gz out=${out} chrom=${chr} >impute_chr${chr}.out &
done

