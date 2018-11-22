#!/bin/bash
# Given two binary PLINK files, merges them
# Syntax: bash merge_bplink_files.sh <dataset1_filename> <dataset2_filename> <output_filename_starting_num_prefix>
# Example: bash merge_bplink_files.sh 03_spit_for_sci_ctrls 38_ceu_cts_IIDupdated_18nov2014update 4
# All formats should be PLINK binary format

prefix=$3

echo "1. Listing Dataset1 SNPs"
cut -f2 $1.bim >temp.dataset_snplist.txt

echo "2. Extracting Dataset common SNPs from dataset2"
plink --noweb --bfile $2 --extract temp.dataset_snplist.txt  --make-bed --out ${prefix}_dataset2_common_snps

echo "3. Listing Dataset2 common SNPs"
cut -f2 ${prefix}_dataset2_common_snps.bim >temp.dataset2_common_snps.txt
prefix_plus1=$[$prefix + 1]

echo "4. Extracting common SNPs from Dataset1"
plink --noweb --bfile $1 --extract temp.dataset2_common_snps.txt --make-bed --out ${prefix_plus1}_dataset1_common_snps
prefix_plus2=$[$prefix_plus1 + 1]

echo "5. Attempting 1st merge"
plink --noweb --bfile ${prefix}_dataset2_common_snps --bmerge ${prefix_plus1}_dataset1_common_snps.bed ${prefix_plus1}_dataset1_common_snps.bim ${prefix_plus1}_dataset1_common_snps.fam --make-bed --out ${prefix_plus2}_mergedsets1
prefix_plus3=$[$prefix_plus2 + 1]

echo "6. Flip strand at specific identified SNPs"
plink --noweb --bfile ${prefix_plus1}_dataset1_common_snps --flip ${prefix_plus2}_mergedsets1.missnp --make-bed --out ${prefix_plus3}_dataset1_flipped_snps
prefix_plus4=$[$prefix_plus3 + 1]

echo "7. Attempting 2nd final merge"
plink --bfile ${prefix_plus3}_dataset1_flipped_snps --noweb --bmerge ${prefix}_dataset2_common_snps.bed ${prefix}_dataset2_common_snps.bim ${prefix}_dataset2_common_snps.fam --make-bed --out ${prefix_plus3}_mergedsets2

echo "Done"
