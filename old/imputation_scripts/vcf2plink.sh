#!/bin/bash
# Script to convert a BEAGLE-imputed VCF file and convert it to PLINK 
# By default, takes all SNPs with AR2>0.3
# Syntax: bash vcf2plink.sh <vcf_filename w/out .vcf.gz> <output_filename_prefix>
# Example: Say we wanted to convert 10_UKRE_chr11.imputed.vcf.gz to PLINK format
#          bash vcf2plink.sh 10_UKRE_chr11.imputed 11_UKRE_chr11.imputed.AR2.0.3

dir="/home/naim/UKRE_imputation/150317-561UKRE_CEU_FULL/imputation_output/"
#minAR2=0.3

echo "1. Extracting AR2"
zcat ${dir}$1.vcf.gz |grep -v ^# |cut -f8 |cut -d';' -f1 |cut -d'=' -f2 >temp.AR2.$1
echo "2. Removing VCF header"
zcat ${dir}$1.vcf.gz |grep -v ^# >temp.$1
echo "3. Adding AR2 field"
paste temp.AR2.$1 temp.$1 >temp2.$1

echo "Some cleanup"
rm temp.AR2.$1
rm temp.$1

echo "4. Filtering out SNPs below 0.3"
awk '$1 > 0.7' temp2.$1 |cut -f2- >temp3.$1

echo "Some cleanup"
rm temp2.$1

echo "5. Extracting original header"
zcat ${dir}$1.vcf.gz |grep ^# >temp.header.$1
echo "6. Re-assembling into new VCF file"
cat temp.header.$1 temp3.$1 |bgzip >temp4.$1.vcf.gz

echo "Some cleanup"
rm temp3.$1
rm temp.header.$1

echo "7. Removing indels using vcftools"
vcftools --gzvcf temp4.$1.vcf.gz --remove-indels --recode --recode-INFO-all --out temp4.snps_only.$1

echo "Some cleanup"
rm temp4.$1.vcf.gz

echo "8. Converting to PLINK format using vcftools"
vcftools --vcf temp4.snps_only.$1.recode.vcf --plink --out $2

echo "Some cleanup"
rm temp4.snps_only.$1.recode.vcf



#echo "9. Correcting FID column in PLINK ped file (this is only specific to Rolandic Epilepsy FID/IID coding)"
#echo "9.1 Generating correct FIDs"
#awk '{gsub("_...","",$1); print $1}' $2.ped >temp.fids.$2.txt
#echo "9.2 Cutting out FIDs"
#cut -f2- $2.ped >temp.cutped.$2.txt
#echo "9.3 Making it space-delimited rather than tab-delimited"
#tr "\t" " " <temp.cutped.$2.txt >temp.cutped2.$2.txt
#echo "9.4 Pasting it all together"
#paste -d' ' temp.fids.$2.txt temp.cutped2.$2.txt >$2.ped
