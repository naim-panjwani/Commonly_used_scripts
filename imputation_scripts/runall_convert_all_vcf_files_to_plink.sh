#!/bin/bash
# This script will simply send out the jobs for conversion of imputed VCF files to PLINK format
# Syntax: bash convert_all_vcf_files_to_plink.sh <output_filename_prefix1> <output_filename_prefix2>
# Example: bash convert_all_vcf_files_to_plink.sh UKRE.exome
# Last modified: 31-Mar-2013

OUTVAR=$1
OUTVAR2=$2

for chri in {1..23}; do
  chr=$chri
  if [[ $chri == "23" ]]; then chr="X"; fi
    echo "Starting conversion of Chromosome ${chr}"
    jobid=`qsub -l mem=128g,vmem=128g,nodes=1:ppn=1,walltime=72:00:00 \
                -o /home/naim/UKRE_imputation/150317-561UKRE_CEU_FULL/imputation_output/vcf2plink/joboutdir/ \
                -e /home/naim/UKRE_imputation/150317-561UKRE_CEU_FULL/imputation_output/vcf2plink/joboutdir/ \
                -d /home/naim/UKRE_imputation/150317-561UKRE_CEU_FULL/imputation_output/vcf2plink/ \
                -N chr${chr} \
                -v chr=${chr},OUTVAR=${OUTVAR},OUTVAR2=${OUTVAR2} \
                   convert_all_vcf_files_to_plink.pbs`
done
