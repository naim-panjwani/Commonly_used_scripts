#!/bin/bash
# Takes a VCF file whose CHROM field precedes by "chr" and removes "chr" prefix, plus replaces chromosome "23" labels with "X"
# Also splits the original file into one file per chromosome
# This is in order to make it compatible with the 1000 Genomes reference files to use with BEAGLE
# Syntax: bash repair_chr_field.sh <filename w/out the .vcf.gz extension> <outfile w/out .vcf.gz extension>
# Example: bash repair_chr_field.sh 03_UKRE_Spit_for_Sci_unwanted_snps_rm 04_UKRE_Spit_for_Sci_unwanted_snps_rm

sample=$1
samplename="${sample}.vcf.gz"
outname=$2

for i in {1..22}; do
  echo "Chr ${i}"
  outfile="${outname}.chr${i}.vcf.gz"
  (zcat ${samplename} |grep ^"#"; zcat ${samplename} |grep -e ^chr${i}\\s |sed 's/^chr//g') |bgzip >${outfile}
done
echo "Chr X"
outfile="${outname}.chrX.vcf.gz"
(zcat ${samplename} |grep ^"#"; zcat ${samplename} |grep ^23 |sed 's/^23/X/g') |bgzip >${outfile}

echo "Done"
