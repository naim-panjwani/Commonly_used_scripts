#!/bin/bash
# This script runs all 5 steps GStoPLINK sequentially
# Dependencies: create_cM_file.py
# Syntax: bash GStoPLINK_runall.sh <FinalReport.txt file> <SNPTable.txt file> <snp column in SNPTable_no_header.txt> <chromosome column> <snp coordinates column> <Sample_Table.txt filename> <number of samples> <IID column> <gender column> 

bash GStoPLINK_step1.sh $1 $2
bash GStoPLINK_step2.sh
bash GStoPLINK_step3.sh $3 $4 $5
bash GStoPLINK_step4.sh
bash GStoPLINK_step5.sh $5 $6 $7 $8

echo "Intermediary files created (list stored under intermediary_files_list.txt):"
echo "FinalReportHeader.txt
FinalReport_no_header.txt
SNPTable_no_header.txt
FinalReport_snp_order.txt
SNPTable_snp_order.txt
Step3_cM
Step3_snps
Step3_chr
Step3_loc
Step3_mapfile
genotype1.txt
genotype2.txt
genotype3.txt
iids1_finalreport.txt
iids2_finalreport.txt
iids_sampletable.txt
fids.txt
gender.txt
pheno.txt
pid.txt"  >intermediary_files_list.txt

