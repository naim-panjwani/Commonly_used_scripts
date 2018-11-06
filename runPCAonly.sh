#!/bin/bash
# Prepares parfile with outlier removal and runs PCA 
# Syntax: bash runPCAonly.sh <pruned_dataset_binary_plink> <pruned_dataset_linkage_plink> <output_filename_prefix_num>
# Example: bash runPCAonly.sh 15_pruned_mergedset_1iter_outliers_rm 16_pruned_mergedset_1iter_outliers_rm 17 

prefix=$3
if [[ $prefix < 10 ]]; then
  char_prefix="0${prefix}"
else
  char_prefix=$prefix
fi

echo "1. Creating ${char_prefix}_parfile-no-outlier-removal.txt"
echo "genotypename: $2.ped
snpname: $2.map
indivname: $1.fam
evecoutname: ${char_prefix}_dataset_outlier_rm.eigenvector
evaloutname: ${char_prefix}_dataset_outlier_rm.eigenvalue
numoutlieriter: 0
snpweightoutname: ${char_prefix}_dataset_outlier_rm.snp.eigenvector
hecksizemode: NO" >${char_prefix}_parfile-no-outlier-removal.txt

echo "2. Starting PCA"
#read cont
#if [[ $cont = "y" ]]; then
  nohup smartpca -p ${char_prefix}_parfile-no-outlier-removal.txt >${char_prefix}_pca.out &
  echo "SmartPCA started"
#else
#  echo "Not started"
#fi

