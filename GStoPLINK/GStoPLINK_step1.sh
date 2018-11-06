#!/bin/bash
# $1 Final_Report.txt from GS in Matrix format without the GenCall Score
# $2 SNP_Table.txt from GS
# Example: bash GStoPLINK_step1.sh TEST_Familial_Epilepsy_FinalReport2.txt TEST_Familial_Epilepsy_SNP_Table.txt 

FinalReport=$1
SNPTable=$2

echo "Step 1. Remove header from ${FinalReport} and ${SNPTable}"
head -n 10 $FinalReport >FinalReportHeader.txt
sed '1,10d' $FinalReport >FinalReport_no_header.txt
sed '1d' ${SNPTable} >SNPTable_no_header.txt
echo "Done"
num_snps_FinalReport=`wc -l FinalReport_no_header.txt`
num_snps_SNPTable=`wc -l SNPTable_no_header.txt`
echo "Number of SNPs FinalReport: ${num_snps_FinalReport}"
echo "Number of SNPs SNPTable: ${num_snps_SNPTable}"
echo "Done"

