#!/bin/bash
# Takes pruned PLINK dataset input for KING pca
# Dependency: king_pca_analysis_v3.R
# Syntax: bash runPCA_v5_KING.sh <dataset> <filename_prefix_num>
# Example: bash runPCA_v5_KING.sh 35_merged_hapmap 38

dataset=$1
prefix=$2

echo "Running command:"
echo "bash runPCA_v8_KING.sh ${dataset} {"

echo "1. Listing X chromosome SNPs, indels and ambiguous SNPs to delete"
grep ^23 $1.bim |cut -f2 >snps_to_rm.txt
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' $1.bim |cut -f2 >>snps_to_rm.txt
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' $1.bim |cut -f2 >>snps_to_rm.txt

plink --noweb --bfile $1 --exclude snps_to_rm.txt --make-bed --out ${prefix}_bad_snps_rm


echo "2. Remove long LD region, rare alleles, and prune it"
plink --noweb --bfile ${prefix}_bad_snps_rm --exclude extended_ld_regions.txt --range --make-bed --out ${prefix}_longregionRemoved
#prune
plink --noweb --bfile ${prefix}_longregionRemoved --maf 0.05 --make-bed --out ${prefix}_MAFlongregionRemoved
plink1.9 --noweb --bfile ${prefix}_MAFlongregionRemoved --indep-pairwise 1500 100 0.2 --out ${name}_pruned
plink --noweb --bfile ${prefix}_MAFlongregionRemoved --extract ${name}_pruned.prune.in --make-bed --out ${name}_pruned



echo "3. Starting KING PCA"
king -b ${name}_pruned.bed --pca 20 --prefix ${prefix}_mergedset_pca >${prefix}_king_pca.out
echo "KING PCA DONE"
echo "Generating Eigenvalues file"
grep 'eigenvalues' ${prefix}_king_pca.out |cut -d' ' -f4- |tr ' ' '\n' >${prefix}_eigenvalues.txt
echo "Tracy Widom Analysis"
bash ~/scripts/tracy-widom.sh ${prefix}_eigenvalues.txt ${prefix}_tracy-widom.txt
echo "Replacing X's in PCA results file by NA's"
sed 's/X/NA/g' ${prefix}_mergedset_pcapc.ped >${prefix}tmp
cat ${prefix}tmp >${prefix}_mergedset_pcapc.ped
rm ${prefix}tmp

echo "3. Running a general KING PCA analysis"
Rscript king_pca_analysis_v3.R ${prefix}_eigenvalues.txt ${prefix}_mergedset_pcapc.ped $1.fam ${prefix}_king_analysis 

echo "Script Done"
