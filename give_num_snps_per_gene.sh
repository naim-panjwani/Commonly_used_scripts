#!/bin/bash
# Input: a space-delimited gene list with chr,start,end,gene_name per line and the bim file to count from (no header allowed)
#        the bimfile to count the SNPs from
#        desired output filename 
# Output: table with the number of SNPs per gene
# Syntax: bash give_num_snps_per_gene.sh [gene_coord_filename.txt] [bim_filename.bim] [output_filename.txt]

genes=$1
bimfile=$2
output_name=$3
output_csv="`echo $output_name |cut -d'.' -f1`.csv"
touch $output_name
while read i; do
  chr=`echo $i |cut -d' ' -f1`
  start=`echo $i |cut -d' ' -f2`
  end=`echo $i |cut -d' ' -f3`
  genename=`echo $i |cut -d' ' -f4`
  Rscript count_num_snps.R $chr $start $end $bimfile
  num_snps=`cat tmp_num_snps.txt`
  echo "${chr} ${start} ${end} ${genename} ${num_snps}" >>${output_name}
  echo "${chr},${start},${end},${genename},${num_snps}" >>${output_csv}
done < $genes

