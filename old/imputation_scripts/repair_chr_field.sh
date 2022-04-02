#!/bin/bash
# Takes a VCF file whose CHROM field precedes by "chr" and removes "chr" prefix, plus replaces chromosome "23" labels with "X"
# This is in order to make it compatible with the 1000 Genomes reference files to use with BEAGLE
# Syntax: bash repair_chr_field.sh <filename w/out the .vcf.gz extension>
# Example: bash repair_chr_field.sh 03_UKRE_Spit_for_Sci_unwanted_snps_rm

sample=$1
samplename="${sample}.vcf.gz"
echo "1. Saving header to sample_header.vcf"
gunzip -c ${samplename} |grep '^#' >sample_header.vcf
echo "2. Unzipping and creating file with no header"
gunzip -c ${samplename} |sed '1,4d' >temp.${sample}_noheader.vcf
echo "3. Removing CHROM field"
cut -f2- temp.${sample}_noheader.vcf >temp.${sample}_noheader_nochrom.vcf
echo "4. Separating out CHROM field to separate file"
cut -f1 temp.${sample}_noheader.vcf >temp.${sample}_chromcol.vcf
echo "5. Replacing 23's with X and removing 'chr' prefix"
sed 's/23/X/g' temp.${sample}_chromcol.vcf >temp.${sample}_chromcolX.vcf
sed 's/chr//g' temp.${sample}_chromcolX.vcf >temp.${sample}_chromcolX2.vcf
echo "6. Pasting it all together"
paste temp.${sample}_chromcolX2.vcf temp.${sample}_noheader_nochrom.vcf >temp.${sample}.mod.vcf
echo "7. Putting the header back on and re-zipping"
cat sample_header.vcf temp.${sample}.mod.vcf |bgzip >final.${sample}.vcf.gz
echo "Done"
echo "Check final file final.${sample}.vcf.gz then rename/overwrite to ${sample}.vcf.gz if it worked"
