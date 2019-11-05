#!/bin/bash
# Checks whether the snp orders in FinalReport_no_header.txt and SNPTable_no_header.txt files match
# MAKE SURE SNPs IN FinalReport_no_header.txt ARE IN COLUMN 1 and SNPs in SNPTable_no_header.txt ARE IN COLUMN 2 
# Example: GStoPLINK_step2.sh

FinalReport="FinalReport_no_header.txt"
SNPTable="SNPTable_no_header.txt"

echo "Step2: Checking whether the order of the marker in ${FinalReport} and ${SNPTable} are the same or not"
cut -f1 ${FinalReport} >FinalReport_snp_order.txt
cut -f2 ${SNPTable} >SNPTable_snp_order.txt
diff FinalReport_snp_order.txt SNPTable_snp_order.txt
echo "Done"

