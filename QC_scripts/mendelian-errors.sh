#!/bin/bash

# PLEASE NOTE THAT THERE IS A MANUAL STEP HERE, SO DON'T RUN THIS SCRIPT BLINDLY...


bplink=$1
outprefix=$2

tmp=$(date '+%s')

echo 
echo "Changing phenotype column of -9 to 2 for PEDSTATS compatibility"
echo
cp ${bplink}.fam ${tmp}-${bplink}.fam # creating backup of fam file
sed -i 's/-9/2/g' ${bplink}.fam

echo 
echo "Extracting the families from PLINK file"
echo families=$(cut -d' ' -f1 ${bplink}.fam |sort |uniq -d)
(for family in $families; do grep ^"${family}\s" ${bplink}.fam; done) >${tmp}-keepfile.fam
plink --bfile ${bplink} --keep ${tmp}-keepfile.fam --make-bed --out ${tmp}-fams_only


echo 
echo "Identifying any dummy parents and adding them"
echo 
# MANUAL STEP... Example follows:
echo "Creating dummy ped file"
echo "12 Fake13 0 0 2 1 A A" >${tmp}-dummy_parents.ped
echo "21 Fake11 0 0 1 1 A A" >>${tmp}-dummy_parents.ped
echo "21 Fake12 0 0 2 1 A A" >>${tmp}-dummy_parents.ped
echo "Creating dummy map file"
echo "1 rsnp 0 1000000" >${tmp}-dummy_parents.map


echo 
echo "Merging dummy parents manually created"
echo 
plink --bfile ${tmp}-fams_only --merge ${tmp}-dummy_parents --make-bed --out ${tmp}-dummy_parents_merge

echo 
echo "Extracting autosomal SNPs and chrX SNPs into separate files"
echo
dataset=${tmp}-dummy_parents_merge 
awk '$1 < 1 || $1 > 22' ${dataset}.bim |cut -f2 >${tmp}-snps_to_rm.txt
grep ^23 ${dataset}.bim |cut -f2 >${tmp}-chrX_snps.txt
cat ${tmp}-chrX_snps.txt >>${tmp}-snps_to_rm.txt
plink --bfile ${dataset} --exclude ${tmp}-snps_to_rm.txt --recode --out ${tmp}-autosome
plink --bfile ${dataset} --extract ${tmp}-chrX_snps.txt --recode --out ${tmp}-chrX


echo 
echo "Creating dat files for autosome and X chromosome"
echo 
(echo -e "A\tB"; awk '{print "M\t"$2}' ${tmp}-autosome.map) >${tmp}-autosome.dat
(echo -e "A\tB"; awk '{print "M\t"$2}' ${tmp}-chrX.map) >${tmp}-chrX.dat

echo 
echo "Running PEDSTATS"
echo 

tmp2=$(( $tmp + 1 ))
pedstats -d ${tmp}-autosome.dat -p ${tmp}-autosome.ped >${tmp2}-pedstats.log
pedstats -d ${tmp}-chrX.dat -p ${tmp2}-chrX.ped --chromosomeX >${tmp2}-chrX_pedstats.log

grep "\- Fam" ${tmp2}-pedstats.log |cut -d' ' -f1 >${tmp2}-mendel_error_snps_autosomes.txt
grep "\- Fam" ${tmp2}-chrX_pedstats.log |cut -d' ' -f1 >${tmp2}-mendel_error_snps_chrX.txt

cat ${tmp2}-mendel_error_snps_autosomes.txt ${tmp2}-mendel_error_snps_chrX.txt > ${outprefix}-Mendel_error_SNPs.txt

rm ${tmp}*
rm ${tmp2}*



