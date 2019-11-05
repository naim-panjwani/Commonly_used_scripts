#!/bin/bash
# Given the dataset and the hapmap reference, merges them
# Syntax: bash merge_with_hapmap.sh <dataset_filename> <hapmap_filename> <output_filename_starting_num_prefix>
# Example: bash merge_with_hapmap.sh Step31_CTS_cases_23relateds 07_hapmap_hg19_related_rm_ceu_tsi 32
# Hapmap data file should be under ./hapmap_data
# All formats should be PLINK binary format

hapmap_dir="hapmap_data/"
prefix=$3

echo "1. Listing Dataset SNPs"
cut -f2 $1.bim >temp.dataset_snplist.txt

echo "2. Extracting Dataset common SNPs from Hapmap reference"
plink --noweb --bfile ${hapmap_dir}$2 --extract temp.dataset_snplist.txt  --make-bed --out ${prefix}_hapmap_common_snps

echo "3. Listing Hapmap common SNPs"
cut -f2 ${prefix}_hapmap_common_snps.bim >temp.hapmap_common_snps.txt
prefix_plus1=$[$prefix + 1]

echo "4. Extracting common SNPs from Dataset"
plink --noweb --bfile $1 --extract temp.hapmap_common_snps.txt --make-bed --out ${prefix_plus1}_dataset_hapmap_common_snps
prefix_plus2=$[$prefix_plus1 + 1]

echo "5. Attempting 1st merge"
plink --noweb --bfile ${prefix}_hapmap_common_snps --bmerge ${prefix_plus1}_dataset_hapmap_common_snps.bed ${prefix_plus1}_dataset_hapmap_common_snps.bim ${prefix_plus1}_dataset_hapmap_common_snps.fam --make-bed --out ${prefix_plus2}_merged_hapmap
prefix_plus3=$[$prefix_plus2 + 1]

echo "6. Flip strand at specific identified SNPs"
plink --noweb --bfile ${prefix_plus1}_dataset_hapmap_common_snps --flip ${prefix_plus2}_merged_hapmap.missnp --make-bed --out ${prefix_plus3}_dataset_flipped_snps
prefix_plus4=$[$prefix_plus3 + 1]

echo "7. Attempting 2nd final merge"
plink --bfile ${prefix_plus3}_dataset_flipped_snps --noweb --bmerge ${prefix}_hapmap_common_snps.bed ${prefix}_hapmap_common_snps.bim ${prefix}_hapmap_common_snps.fam --make-bed --out ${prefix_plus3}_merged_hapmap

echo "Done"
 
