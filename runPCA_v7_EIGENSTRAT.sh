#!/bin/bash
# Takes pruned PLINK dataset input for EIGENSTRAT PCA
# Syntax: bash runPCA_v7_EIGENSTRAT.sh <pruned_dataset> <filename_prefix_num>
# Example: bash runPCA_v7_EIGENSTRAT.sh 35_merged_hapmap_pruned 38

prefix=$2

echo "1. Listing X chromosome SNPs, indels and ambiguous SNPs to delete"
grep ^23 $1.bim |cut -f2 >snps_to_rm.txt
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' $1.bim |cut -f2 >>snps_to_rm.txt
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' $1.bim |cut -f2 >>snps_to_rm.txt
num_ind=`wc -l $1.fam |cut -d' ' -f1`
awk '{print $1,$2}' $1.fam >tmp1${prefix}
yes "1" |head -n $num_ind >tmp2${prefix}
paste -d' ' tmp1${prefix} tmp2${prefix} >${prefix}_pheno.txt
rm tmp1${prefix} tmp2${prefix}

plink --noweb --bfile $1 --exclude snps_to_rm.txt --pheno ${prefix}_pheno.txt --make-bed --out ${prefix}_bad_snps_rm
plink --noweb --bfile ${prefix}_bad_snps_rm --recode --out ${prefix}_bad_snps_rm_recode

echo "2. Generating parfile"
echo "genotypename: ${prefix}_bad_snps_rm_recode.ped
snpname: ${prefix}_bad_snps_rm_recode.map
indivname: ${prefix}_bad_snps_rm.fam
evecoutname: ${prefix}_eigenpcapc.eigenvector
evaloutname: ${prefix}_eigenpcapc.eigenvalue
numoutlieriter: 0
snpweightoutname: ${prefix}_eigenpcapc_snpweights.eigenvector
hecksizemode: NO" >${prefix}_smartpca_parfile.txt

echo "3. Starting EIGENSTRAT PCA"
smartpca -p ${prefix}_smartpca_parfile.txt >${prefix}_smartpca.out
echo "EIGENSTRAT PCA DONE"
echo "Tracy Widom Analysis"
bash ~/scripts/tracy-widom.sh ${prefix}_eigenpcapc.eigenvalue ${prefix}_tracy-widom.txt

#echo "3. Running a general EIGENSTRAT PCA analysis"
#Rscript king_pca_analysis_v3.R ${prefix}_eigenvalues.txt ${prefix}_mergedset_pcapc.ped $1.fam ${prefix}_king_analysis 

echo "Script Done"
