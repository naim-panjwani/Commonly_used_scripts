#!/bin/bash
# Script to convert a BEAGLE-imputed VCF file and convert it to PLINK 
# By default, takes all SNPs with AR2>0.3
# Syntax: bash vcf2plink.sh <vcf_filename w/out .vcf.gz> <output_filename>
# Example: Say we wanted to convert 10_UKRE_chr11.imputed.vcf.gz to PLINK format
#          bash vcf2plink.sh 10_UKRE_chr11.imputed 11_UKRE_chr11.imputed.AR2.0.3
 
set -e
TEMPDIR="/localhd/${PBS_JOBID}"

echo "1. Extracting AR2"
zcat $1.vcf.gz |grep -v ^# |cut -f8 |cut -d';' -f1 |cut -d'=' -f2 >$TEMPDIR/temp.AR2.$1
echo "2. Removing VCF header"
zcat $1.vcf.gz |grep -v ^# >$TEMPDIR/temp.$1
echo "3. Adding AR2 field"
paste $TEMPDIR/temp.AR2.$1 $TEMPDIR/temp.$1 >$TEMPDIR/temp2.$1

echo "Some cleanup"
rm $TEMPDIR/temp.AR2.$1
rm $TEMPDIR/temp.$1

echo "4. Filtering out SNPs below 0.8"
awk '$1 > 0.8' $TEMPDIR/temp2.$1 |cut -f2- >$TEMPDIR/temp3.$1

echo "Some cleanup"
rm $TEMPDIR/temp2.$1

echo "5. Extracting original header"
zcat $1.vcf.gz |grep ^# >$TEMPDIR/temp.header.$1
echo "6. Re-assembling into new VCF file"
cat $TEMPDIR/temp.header.$1 $TEMPDIR/temp3.$1 |bgzip >$TEMPDIR/temp4.$1.vcf.gz

echo "Some cleanup"
rm $TEMPDIR/temp3.$1
rm $TEMPDIR/temp.header.$1

echo "7. Recoding using VCFTools"
vcftools --gzvcf $TEMPDIR/temp4.$1.vcf.gz --recode --recode-INFO-all --out $TEMPDIR/temp4.snps_only.$1

echo "Some cleanup"
rm $TEMPDIR/temp4.$1.vcf.gz

echo "8. Converting to PLINK format using vcftools"
vcftools --vcf $TEMPDIR/temp4.snps_only.$1.recode.vcf --plink --out $TEMPDIR/temp.$2

echo "Some cleanup"
rm $TEMPDIR/temp4.snps_only.$1.recode.vcf

echo "9. Converting to binary PLINK format"
plink --noweb --file $TEMPDIR/temp.$2 --make-bed --out $2

echo "Some cleanup"
rm $TEMPDIR/temp.$2*

echo "Done"

