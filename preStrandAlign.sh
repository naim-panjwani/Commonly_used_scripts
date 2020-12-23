#!/bin/bash
# PRE: Takes PLINK binary fileset
# POST: Removes monomorphic SNPs, indels, and 
#       SNPs not in chromosomes 1-23 (keeps chr25 pseudoautosomal region of X also)
#       Also converts underscores in FID/IID to dashes
# Example: bash preStrandAligh.sh 41-JME_dup_snps_rm 44-JME_unwanted_snps_rm 44
# NOTE: You may also want to remove unmappable SNPs prior to strand alignment (listed .miss and .multiple Will Rayner files)

plinkfile=$1 # enter without the .bed/.bim/.fam extensions
outfilename=$2
prefix=$3

# Remove monomorphic SNPs first:
echo "Gathering monomorphic SNPs"
awk '$5==0' ${plinkfile}.bim >${prefix}_snps_to_rm.txt
num_snps=$(wc -l ${prefix}_snps_to_rm.txt |cut -d' ' -f1)
echo "There are " ${num_snps} " monomorphic SNPs"
echo "Gathering indels (I/D SNPs)"
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' ${plinkfile}.bim >${prefix}_temporary
cat ${prefix}_temporary >>${prefix}_snps_to_rm.txt
num_snps=$(wc -l ${prefix}_temporary |cut -d' ' -f1)
echo "There are " ${num_snps} " indels"
rm ${prefix}_temporary
echo -e "Gathering SNPs not in chromosomes 1-23;\nWe also want \"chr25\" which plink codes for the pseudoautosomal region of X"
awk '$1 < 1 || $1 > 25 || $1 == 24' ${plinkfile}.bim >${prefix}_temporary
cat ${prefix}_temporary >>${prefix}_snps_to_rm.txt
num_snps=$(wc -l ${prefix}_temporary |cut -d' ' -f1)
echo "There are " ${num_snps} " non-chr1-23 SNPs"
rm ${prefix}_temporary

echo 
echo "Removing identified SNPs with PLINK, and changing underscores to dashes in ID fields"
echo 
# Underscores in FID/IID fields will cause grief later on, so change it to a dash:
paste -d' ' <(cut -d' ' -f1,2 ${plinkfile}.fam) <(cut -d' ' -f1,2 ${plinkfile}.fam |sed 's/_/-/g') >${prefix}_ID_underscore_to_dash.txt
plink --bfile ${plinkfile} --exclude ${prefix}_snps_to_rm.txt --update-ids ${prefix}_ID_underscore_to_dash.txt --make-bed --out ${outfilename}

echo -e "Done\n"
echo "NOTE: You may also want to remove unmappable SNPs prior to strand alignment (listed .miss and .multiple Will Rayner files)"

