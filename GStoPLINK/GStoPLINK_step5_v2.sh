#!/bin/bash
# Final step to create tfam file
# Syntax: bash GStoPLINK_step5.sh <Sample_Table.txt filename> <number of samples> <IID column> <gender column>
# Example for 288 individuals, IIDs on column 2 of SampleTable, and genders in column 6 
# bash GStoPLINK_step5.sh TEST_Familial_Epilepsy_Sample_Table.txt 288 2 6


FinalReport="FinalReport_no_header.txt"
SNPTable="SNPTable_no_header.txt"
FinalReportHeader="FinalReportHeader.txt"
SamplesTable=$1

echo "Step 5. Generating tfam file"
#echo "Step 5.1 Generating iids in FinalReport.txt"
#head -1 ${FinalReportHeader} >iids1_finalreport.txt
#cat iids1_finalreport.txt |tr "\t" "\n" |sed '1d' |tr "-" "_" |sed 's/\r//g' >iids2_finalreport.txt
#echo "Step 5.2 Generating iids in SampleTable"
#sed '1d' ${SamplesTable} |cut -f$3 |tr "-" "_" >iids_sampletable.txt
#echo "Step 5.3 Comparing IID order in ${FinalReport} and ${SamplesTable}"
#diff iids2_finalreport.txt iids_sampletable.txt
##cat iids2_finalreport.txt |tr "_" "\t" |cut -f1 >fids.txt
#echo "Step 5.4 Generating fids.txt"
#cat iids2_finalreport.txt >fids.txt
#echo "Step 5.5 Generating gender.txt"
#sed '1d' ${SamplesTable} |cut -f$4 |sed 's/M/1/g' |sed 's/F/2/g' >gender.txt
#numsamples=`wc -l iids2_finalreport.txt` |cut -d' ' -f1
#echo "Step 5.6 Creating pid.txt file with $2 zero lines and pheno.txt file with -9 lines"
#python create_cM_file.py -i $2 -o pid.txt
#python create_pheno_file.py -i $2 -o pheno.txt
#echo "Step 5.7 Generating pre-final genodata.tfam file"
#paste -d' ' fids.txt iids_sampletable.txt pid.txt pid.txt gender.txt >genodata.tfam
#echo "Step 5.8 Adding -9 phenotype column using R"
#R CMD BATCH add_pheno_column.R
#echo "Created genodata.tfam file"

sed '1d' $SamplesTable |awk '{print $1"\t"$1"\t"0"\t"0"\t"$2"\t"0}' >genodata.tfam

