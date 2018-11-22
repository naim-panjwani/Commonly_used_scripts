#!/bin/bash
# Take a binary PLINK file and splits it into the 23 chromosomes; destination files will be .raw PLINK files (recodeA)
# Syntax: bash split_by_chrom.sh <source filename> <destination filename> 
# Example: bash split_by_chrom.sh 15_mergedset_pca_outliers_rm 20_mergedset_exome_cts_and_spit_for_science

outname=$2
for i in {1..23}; do
  nohup plink --noweb --bfile $1 --chr ${i} --make-bed --out ${outname}_bchr${i} &
  nohup plink --noweb --bfile $1 --chr ${i} --recodeA --out ${outname}_chr${i} &
done

