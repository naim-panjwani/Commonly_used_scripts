#!/bin/bash
# Will simply send out the jobs for conversion of imputed VCF files to PLINK format
# Syntax: bash convert_all_vcf_files_to_plink.sh <output_filename_prefix>
# Example: bash convert_all_vcf_files_to_plink.sh UKRE.exome

for i in {1..22}
  do
    echo "Starting conversion of Chromosome ${i}"
    bash vcf2plink.sh UKRE_data.chr${i}.imputed $1.chr${i}.imputed.AR2.0.3 >chr${i}.out
    echo "Removing temporary files"
    rm temp*
  done
echo "Starting conversion of Chromosome X"
bash vcf2plink.sh UKRE_data.chrX.imputed $1.chrX.imputed.AR2.0.3 >chrX.out
echo "Removing temporary files"
rm temp*
