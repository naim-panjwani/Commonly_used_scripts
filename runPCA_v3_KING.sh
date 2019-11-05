#!/bin/bash
# Prepares pruned PLINK dataset input for smartpca and then runs PCA
# Syntax: bash runPCA_v3_KING.sh <unpruned_dataset> <prune.in_snps> <filename_prefix_num>
# Example: bash runPCA_v3_KING.sh 35_merged_hapmap 36_CTS_relateds_hapmap_pruned 38

prefix=$3
if [[ $prefix < 10 ]]; then
  char_prefix=0${prefix}
else
  char_prefix=$prefix
fi

echo "1. Extracting pruned in SNPs"
plink --noweb --bfile $1 --extract $2.prune.in --make-bed --out ${char_prefix}_pruned

prefix_plus1=$[$prefix + 1]
if [[ $prefix_plus1 < 10 ]]; then
  char_prefix1=0${prefix_plus1}
else
  char_prefix1=$prefix_plus1
fi


echo "2. Removing X chromosome SNPs"
grep ^23 ${char_prefix1}_pruned.bim |cut -f2 >${char_prefix1}_chr23snps.txt
plink --noweb --bfile ${char_prefix}_pruned --exclude ${char_prefix1}_chr23snps.txt --make-bed --out ${char_prefix1}_autosomalSNPs

prefix_plus2=$[$prefix_plus1 + 1]
if [[ $prefix_plus2 < 10 ]]; then
  char_prefix2=0${prefix_plus2}
else
  char_prefix2=$prefix_plus2
fi


echo "3. Deleting indels and ambiguous SNPs"
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' ${char_prefix1}_autosomalSNPs.bim |cut -f2 >${char_prefix2}_DI_snps.txt
# List ambiguous A/T G/C SNPs
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $5 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' ${char_prefix1}_autosomalSNPs.bim |cut -f2 >${char_prefix2}_ambiguous_snps.txt
cat ${char_prefix2}_DI_snps ${char_prefix2}_ambiguous_snps.txt >${char_prefix2}_DI_and_ambiguous_snps_to_delete.txt
plink --noweb --bfile ${char_prefix1}_autosomalSNPs --exclude ${char_prefix2}_DI_and_ambiguous_snps_to_delete.txt --make-bed --out ${char_prefix2}_autosomalSNPs_DI_ambiguous_snps_rm


prefix_plus3=$[$prefix_plus2 + 1]
if [[ $prefix_plus3 < 10 ]]; then
  char_prefix3=0${prefix_plus3}
else
  char_prefix3=$prefix_plus3
fi


#echo "4. Creating map and ped files for smartpca"
#plink --noweb --bfile ${char_prefix2}_dataset_hapmap_autosomalSNPs_DI_rm --recode --out ${char_prefix3}_dataset_hapmap_autosomalSNPs_DI_rm

#prefix_plus4=$[$prefix_plus3 + 1]
#if [[ $prefix_plus4 < 10 ]]; then
#  char_prefix4=0${prefix_plus4}
#else
#  char_prefix4=$prefix_plus4
#fi


#echo "5. Creating ${char_prefix4}_parfile-no-outlier-removal.txt"
#echo "genotypename: ${char_prefix3}_dataset_hapmap_autosomalSNPs_DI_rm.ped
#snpname: ${char_prefix3}_dataset_hapmap_autosomalSNPs_DI_rm.map
#indivname: ${char_prefix2}_dataset_hapmap_autosomalSNPs_DI_rm.fam
#evecoutname: ${char_prefix4}_dataset_hapmap_autosomalSNPs_DI_rm.eigenvector
#evaloutname: ${char_prefix4}_dataset_hapmap_autosomalSNPs_DI_rm.eigenvalue
#numoutlieriter: 0
#snpweightoutname: ${char_prefix4}_dataset_hapmap_autosomalSNPs_DI_rm.snp.eigenvector
#hecksizemode: NO" >${char_prefix4}_parfile-no-outlier-removal.txt

#echo "6. Changing phenotype column to 1 in fam file"
#sed 's/-9/1/g' ${char_prefix2}_dataset_hapmap_autosomalSNPs_DI_rm.fam >temp.fam
#cat temp.fam >${char_prefix2}_dataset_hapmap_autosomalSNPs_DI_rm.fam

#echo "7. Starting PCA"
#read cont
#if [[ $cont = "y" ]]; then
#  nohup smartpca -p ${char_prefix4}_parfile-no-outlier-removal.txt >${char_prefix4}_pca.out &
#  echo "SmartPCA started"
#else
#  echo "Not started"
#fi


echo "4. Starting PCA"
king -b ${char_prefix2}_autosomalSNPs_DI_ambiguous_snps_rm.bed --pca --prefix ${char_prefix3}_mergedset_pca >${char_prefix3}_king_pca.out
echo "SmartPCA started"

