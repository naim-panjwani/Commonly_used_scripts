#!/bin/bash
# Takes pruned PLINK dataset input for KING pca
# Dependency: king_pca_analysis_v3.R
# Syntax: bash runPCA_v5_KING.sh <pruned_dataset> <filename_prefix_num>
# Example: bash runPCA_v5_KING.sh 35_merged_hapmap_pruned 38

prefix=$2

echo "1. Listing X chromosome SNPs, indels and ambiguous SNPs to delete"
grep ^23 $1.bim |cut -f2 >snps_to_rm.txt
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' $1.bim |cut -f2 >>snps_to_rm.txt
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' $1.bim |cut -f2 >>snps_to_rm.txt

plink --noweb --bfile $1 --exclude snps_to_rm.txt --make-bed --out ${prefix}_bad_snps_rm


echo "2. Starting KING PCA"
king -b ${prefix}_bad_snps_rm.bed --pca 10 --prefix ${prefix}_mergedset_pca >${prefix}_king_pca.out
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
