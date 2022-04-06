#!/bin/bash
# Given two binary PLINK files, merges them at the common SNPs only
# This version deletes the intermediary files and outputs the merged "output_filename"
# Syntax: bash merge_bplink_files.sh <dataset1_filename> <dataset2_filename> <output_filename> <temp_identifier>
# Example: bash merge_bplink_files_v2.sh 03_spit_for_sci_ctrls 38_ceu_cts_IIDupdated_18nov2014update 04_UK_SS_merged chr11
# All formats should be PLINK binary format

#set -e

outname=$3
id=$4

echo "1. Listing Dataset1 SNPs"
cut -f2 $1.bim >temp_${id}.dataset_snplist.txt

echo "2. Extracting Dataset common SNPs from dataset2"
plink --noweb --bfile $2 --extract temp_${id}.dataset_snplist.txt  --make-bed --out temp_${id}_dataset2_common_snps --memory 4000

echo "3. Listing Dataset2 common SNPs"
cut -f2 temp_${id}_dataset2_common_snps.bim >temp_${id}.dataset2_common_snps.txt

echo "4. Extracting common SNPs from Dataset1"
plink --noweb --bfile $1 --extract temp_${id}.dataset2_common_snps.txt --make-bed --out temp_${id}_dataset1_common_snps --memory 4000

echo "5. Attempting 1st merge"
plink --noweb --bfile temp_${id}_dataset2_common_snps --bmerge temp_${id}_dataset1_common_snps.bed temp_${id}_dataset1_common_snps.bim temp_${id}_dataset1_common_snps.fam --make-bed --out temp_${id}_mergedsets1 --memory 4000

if [[ -e "temp_${id}_mergedsets1-merge.missnp" ]]; then
  echo "6. Flip strand at specific identified SNPs"
  plink --noweb --bfile temp_${id}_dataset1_common_snps --flip temp_${id}_mergedsets1-merge.missnp --make-bed --out temp_${id}_dataset1_flipped_snps --memory 4000

  echo "7. Attempting 2nd merge"
  plink --noweb --bfile temp_${id}_dataset1_flipped_snps --bmerge temp_${id}_dataset2_common_snps.bed temp_${id}_dataset2_common_snps.bim temp_${id}_dataset2_common_snps.fam --make-bed --out temp_${id}_mergedsets2 --memory 4000
  
  if [[ -e "temp_${id}_mergedsets2-merge.missnp" ]]; then
    mv temp_${id}_mergedsets2-merge.missnp polymorphic_snp_list_${id}.missnp
    echo "Likely have polymorphic SNPs; removing polymorphic SNPs listed in polymorphic_snp_list_${id}.missnp"
    plink --noweb --bfile temp_${id}_dataset1_flipped_snps --exclude polymorphic_snp_list_${id}.missnp --make-bed --out temp_${id}_dataset1_polymor_rm --memory 4000
    plink --noweb --bfile temp_${id}_dataset2_common_snps --exclude polymorphic_snp_list_${id}.missnp --make-bed --out temp_${id}_dataset2_polymor_rm --memory 4000
  
    echo "8. Attempting 3rd merge"
    plink --noweb --bfile temp_${id}_dataset1_polymor_rm --bmerge temp_${id}_dataset2_polymor_rm.bed temp_${id}_dataset2_polymor_rm.bim temp_${id}_dataset2_polymor_rm.fam --make-bed --out temp_${id}_mergedsets3 --memory 4000

    if [[ ! -e "temp_${id}_dataset2_polymor_rm.bed" ]]; then
      echo "Something went wrong; check temp_${id}_mergedsets3.log file"
      exit 3
    else
      mv temp_${id}_mergedsets3.bed ${outname}.bed
      mv temp_${id}_mergedsets3.bim ${outname}.bim
      mv temp_${id}_mergedsets3.fam ${outname}.fam
      mv temp_${id}_mergedsets3.log ${outname}.log
      rm temp_${id}*
    fi
  else
    mv temp_${id}_mergedsets2.bed ${outname}.bed
    mv temp_${id}_mergedsets2.bim ${outname}.bim
    mv temp_${id}_mergedsets2.fam ${outname}.fam
    mv temp_${id}_mergedsets2.log ${outname}.log
    rm temp_${id}*
  fi

elif [[ ! -e "temp_${id}_mergedsets1.bed" ]]; then
  echo "Something went wrong; check temp_${id}_mergedsets1.log file"
  echo "Intermediate temporary files have not been deleted"
  exit 1

else
  mv temp_${id}_mergedsets1.bed ${outname}.bed
  mv temp_${id}_mergedsets1.bim ${outname}.bim
  mv temp_${id}_mergedsets1.fam ${outname}.fam
  mv temp_${id}_mergedsets1.log ${outname}.log
  rm temp_${id}*
fi

echo "Done"
