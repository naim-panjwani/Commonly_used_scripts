#!/bin/bash
# Creates the mapfile
# Syntax: bash GStoPLINK_step3.sh <snp column in SNPTable_no_header.txt> <chromosome column> <snp coordinates column>
# Example: bash GStoPLINK_step3.sh 2 3 4

FinalReport="FinalReport_no_header.txt"
SNPTable="SNPTable_no_header.txt"

echo "Step 3: Use SNP file to generate map file"
numsnps=`wc -l ${SNPTable} |cut -d' ' -f1`
echo "There are ${numsnps} SNPs" 

awk '{print $2"\t"$1"\t"0"\t"$3}' $SNPTable >Step3_mapfile
echo "Done"

