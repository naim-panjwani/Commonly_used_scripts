#!/bin/bash
# Prepares pruned PLINK dataset input for smartpca and then runs PCA
# Syntax: bash runPCA.sh <unpruned_dataset> <prune.in_snps> <filename_prefix_num>
# Example: bash runPCA.sh 35_merged_hapmap 36_CTS_relateds_hapmap_pruned 38

prefix=$3

echo "1. Extracting pruned in SNPs"
plink --noweb --bfile $1 --extract $2.prune.in --make-bed --out ${prefix}_dataset_hapmap_pruned

prefix_plus1=$[$prefix + 1]
echo "2. Removing X chromosome SNPs"
grep ^23 ${prefix}_dataset_hapmap_pruned.bim |cut -f2 >${prefix}_chr23snps.txt
plink --noweb --bfile ${prefix}_dataset_hapmap_pruned --exclude ${prefix}_chr23snps.txt --make-bed --out ${prefix_plus1}_dataset_hapmap_autosomalSNPs

prefix_plus2=$[$prefix_plus1 + 1]
echo "3. Deleting indels and ambiguous SNPs"
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' ${prefix_plus1}_dataset_hapmap_autosomalSNPs.bim |cut -f2 >DI_snps.txt
# List ambiguous A/T G/C SNPs
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $5 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' ${prefix_plus1}_dataset_hapmap_autosomalSNPs.bim |cut -f2 >ambiguous_snps.txt
cat DI_snps ambiguous_snps.txt >snps_to_delete.txt
plink --noweb --bfile ${prefix_plus1}_dataset_hapmap_autosomalSNPs --exclude snps_to_delete.txt --make-bed --out ${prefix_plus2}_dataset_hapmap_autosomalSNPs_DI_rm


prefix_plus3=$[$prefix_plus2 + 1]
echo "4. Creating map and ped files for smartpca"
plink --noweb --bfile ${prefix_plus2}_dataset_hapmap_autosomalSNPs_DI_rm --recode --out ${prefix_plus3}_dataset_hapmap_autosomalSNPs_DI_rm

prefix_plus4=$[$prefix_plus3 + 1]
echo "5. Creating ${prefix_plus4}_parfile-no-outlier-removal.txt"
echo "genotypename: ${prefix_plus3}_dataset_hapmap_autosomalSNPs_DI_rm.ped
snpname: ${prefix_plus3}_dataset_hapmap_autosomalSNPs_DI_rm.map
indivname: ${prefix_plus2}_dataset_hapmap_autosomalSNPs_DI_rm.fam
evecoutname: ${prefix_plus4}_dataset_hapmap_autosomalSNPs_DI_rm.eigenvector
evaloutname: ${prefix_plus4}_dataset_hapmap_autosomalSNPs_DI_rm.eigenvalue
numoutlieriter: 0
snpweightoutname: ${prefix_plus4}_dataset_hapmap_autosomalSNPs_DI_rm.snp.eigenvector
hecksizemode: NO" >${prefix_plus4}_parfile-no-outlier-removal.txt

#echo "6. Changing phenotype column to 1 in fam file"
#sed 's/-9/1/g' ${prefix_plus2}_dataset_hapmap_autosomalSNPs_DI_rm.fam >temp.fam
#cat temp.fam >${prefix_plus2}_dataset_hapmap_autosomalSNPs_DI_rm.fam

echo "7. Start PCA?"
read cont
if [[ $cont = "y" ]]; then
  echo "SmartPCA started"
  nohup smartpca -p ${prefix_plus4}_parfile-no-outlier-removal.txt >${prefix_plus4}_pca.out &
else
  echo "Not started"
fi

