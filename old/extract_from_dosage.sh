#!/bin/bash
# Given a compressed dosage file, extracts the specified data between the specified basepairs
# Syntax: bash extract_from_dosage.sh dosage_filename.gz bppositions.txt output_filename
# Example: bash extract_from_dosage.sh chr5/chr5.block1.dosage.gz bppositions.txt temp_test.txt

dosage_dir="/home/data/GWAS2_Consortium/GWAS2_imputation_Oct_2013/GWAS1_imputation/"

bpstarts=`cut -f1 $2`
bpends=`cut -f2 $2`

count=1
touch $3
for i in $bpstarts; do
  current_bpstart=`sed -n "${count}p" $2 |cut -f1`
  current_bpend=`sed -n "${count}p" $2 |cut -f2`
  echo "Extracting ${current_bpstart} to ${current_bpend}"
  zcat ${dosage_dir}$1 |awk '$3>=start && $3<=end' start="$current_bpstart" end="$current_bpend" >>$3 
  count=$[ count + 1 ]
done
