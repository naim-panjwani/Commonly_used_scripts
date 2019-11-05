#!/bin/bash

bplink_file=$1
out_prefix=$2

tmp=$(date '+%s')
echo "Calculating autosomal SNPs heterozygosity"
plink1.9 --bfile ${bplink_file} --het --out ${out_prefix}-autosomes_het
echo "Subsetting chrX SNPs"
plink1.9 --bfile ${bplink_file} --chr 23 --make-bed --out ${out_prefix}-chrX
echo "Hacking chrX to be chr1"
sed 's/^23/1/g' ${out_prefix}-chrX.bim >tmp
cp tmp ${out_prefix}-chrX.bim
echo "Removing temporary files"
rm tmp
echo "Calculating chrX SNPs heterozygosity"
plink1.9 --bfile ${out_prefix}-chrX --het --out ${out_prefix}-chrX_het
echo "Getting the final plot"
Rscript QC_scripts/heterozygosity.R ${out_prefix}-autosomes_het.het ${out_prefix}-chrX_het.het ${bplink_file}.fam ${out_prefix}
