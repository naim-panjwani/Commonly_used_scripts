#!/bin/bash/
# Given the .eigenvector input, removes the header and 1st column to get the twtable
# Syntax: bash make_twtable.sh <.eigenvector_filename> <output_filename_prefix_num>
# Example: bash make_twtable.sh 17_dataset_outlier_rm.eigenvector 18

$prefix=$2
sed '1d' $1 | awk '{out=$2; for(i=3;i<=NF;i++){out=out" "$i}; print out}' >${prefix}_twtable.eigenvector

