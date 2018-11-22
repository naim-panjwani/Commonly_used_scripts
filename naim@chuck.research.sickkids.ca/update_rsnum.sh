#!/bin/bash
# Takes a bim file and a bed-formatted tabix indexed file from dbSNP and outputs a file with
# the subset of rs snp id's corresponding to the base-pair positions in the bim file
# Syntax: bash <bim filename> <dbsnp filename> <output filename>
# Example: bash 38_ceu_cts_IIDupdated_18nov2014update.bim snp141.bed.gz  38_new_bimfile

out=$3
cut -f4 $1 >tempbp_${out}.txt
cut -f1,4 $1 |tr '\t' ':' |paste -d'-' - tempbp_${out}.txt |sed 's/^23:/X:/g' >temppos_${out}.txt
awk '{OFS="" ; print "chr", $1}' temppos_${out}.txt >temppos2_${out}.txt
touch ${out}
for i in `cat temppos2_${out}.txt`
do
  tabix $2 ${i} >tempi_${out}.txt
  echo -e "${i}\tNA\tNA\tNA" >>tempi_${out}.txt #in case no lines are returned
  head -1 tempi_${out}.txt >>${out}
done
rm tempbp_${out}.txt temppos_${out}.txt temppos2_${out}.txt tempi_${out}.txt
