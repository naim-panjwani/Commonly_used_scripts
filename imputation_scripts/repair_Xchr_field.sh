echo "Enter sample file name without the .vcf.gz extension (DON'T FORGET TO ADD .chrX):"
read sample
samplename="${sample}.vcf.gz"
echo "1. Saving header to sample_header.vcf"
gunzip -c ${samplename} |grep '^#' >sample_header.vcf
echo "2. Unzipping and creating file with no header"
gunzip -c ${samplename} |sed '1,4d' >temp.${sample}_noheader.vcf
echo "3. Removing CHROM field"
cut -f2- temp.${sample}_noheader.vcf >temp.${sample}_noheader_nochrom.vcf
echo "4. Separating out CHROM field to separate file"
cut -f1 temp.${sample}_noheader.vcf >temp.${sample}_chromcol.vcf
echo "5. Replacing 23's with X"
sed 's/23/X/g' temp.${sample}_chromcol.vcf >temp.${sample}_chromcolX.vcf
echo "6. Pasting it all together"
paste temp.${sample}_chromcolX.vcf temp.${sample}_noheader_nochrom.vcf >temp.${sample}.mod.vcf
echo "7. Putting the header back on and re-zipping"
cat sample_header.vcf temp.${sample}.mod.vcf |bgzip >final.${sample}.vcf.gz
echo "Done"
echo "Check final file final.${sample}.vcf.gz then rename/overwrite to ${sample}.vcf.gz if it worked"
