#!/bin/bash

plinkfile=$1
outfile=$2

for i in {1..23}; do
  let chrom=$i
  if [[ $i -eq 23 ]]; then chrom="X"; fi
  plink --bfile ${plinkfile} --chr ${i} --recode vcf --out "${outfile}_chr${chrom}" --output-chr MT
  if [[ $i -eq 23 ]]; then # ensuring pseudo-autosomal chrX SNPs make it into the file  
    plink --bfile ${plinkfile} --chr 25 --recode vcf --out ${outfile}_chr25 --output-chr MT
    #sed 's/^23/X/g' ${outfile}_chr${chrom}.vcf >${outfile}_chr${chrom}_fix.vcf
    sed 's/^XY/X/g' ${outfile}_chr25.vcf |grep -v ^"#" >> ${outfile}_chr${chrom}.vcf
    bcftools sort ${outfile}_chr${chrom}.vcf |sed 's/\t0\t/\t0\/0\t/g' | sed 's/\t0\t/\t0\/0\t/g'| sed 's/\t0$/\t0\/0/g'| sed 's/\t1\t/\t1\/1\t/g' | sed 's/\t1\t/\t1\/1\t/g'| sed 's/\t1$/\t1\/1/g' > ${outfile}_chr${chrom}-fix.vcf
    cp ${outfile}_chr${chrom}-fix.vcf ${outfile}_chr${chrom}.vcf
    rm ${outfile}_chr${chrom}-fix.vcf ${outfile}_chr25*
  fi
  bgzip ${outfile}_chr${chrom}.vcf
done

echo "chrX may require further fixing; use fixX.py to fix further"
