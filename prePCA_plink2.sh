#!/bin/bash
# Make sure extended_ld_regions.txt is present locally

dataset=$1
prefix=$2

echo "1. Listing non-autosomal chromosome SNPs, indels and ambiguous SNPs to delete"
awk '$1 < 1 || $1 > 22' ${dataset}.bim |cut -f2 > "${prefix}-snps_to_rm.txt"
#grep ^23 ${dataset}.bim |cut -f2 > "${prefix}-snps_to_rm.txt"
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' ${dataset}.bim |cut -f2 >> "${prefix}-snps_to_rm.txt"
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' ${dataset}.bim |cut -f2 >> "${prefix}-snps_to_rm.txt"
awk '($5 != "A" && $5 != "T" && $5 != "G" && $5 != "C") || ($6 != "A" && $6 != "T" && $6 != "G" && $6 != "C")' ${dataset}.bim |cut -f2 >> "${prefix}-snps_to_rm.txt"

plink2 --bfile ${dataset} --exclude "${prefix}-snps_to_rm.txt" --make-bed --out "${prefix}-bad_snps_rm" --memory 4000


echo "2. Remove long LD region, rare alleles, and prune it"
if [ ! -e extended_ld_regions.txt ]; then ln -s ~/scripts/extended_ld_regions.txt; fi
plink2 --bfile "${prefix}-bad_snps_rm" --exclude range extended_ld_regions.txt --make-bed --out "${prefix}-longregionRemoved" --memory 4000
#prune
plink2 --bfile "${prefix}-longregionRemoved" --maf 0.05 --hwe 0.001 --make-bed --out "${prefix}-MAFlongregionRemoved" --memory 4000
plink2 --bfile "${prefix}-MAFlongregionRemoved" --indep-pairwise 1500 100 0.2 --out "${prefix}-pruned" --memory 4000
plink2 --bfile "${prefix}-MAFlongregionRemoved" --extract "${prefix}-pruned.prune.in" --make-bed --out ${prefix}-pruned --memory 4000

rm "${prefix}-bad_snps_rm"* "${prefix}-longregionRemoved"* "${prefix}-MAFlongregionRemoved"*

