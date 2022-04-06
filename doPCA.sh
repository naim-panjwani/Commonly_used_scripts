#!/bin/bash
# Perform PCA analysis using PC-AiR given a set of individuals and a PLINK file
# Will remove unnecessary SNPs and prune and obtain the kinship matrix via KING prior to running PC-AiR
# Depends on: prePCA_plink2.sh, KING 2.2.4, and PC-AiR_v1_empty_kinFile2.R (tested in conda environment)
# Syntax: bash doPCA.sh <sample_list_to_keep> <plink_root_filename> <outfilename>
# sample_list.txt has FID and IID
# Example: bash ~/scripts/doPCA.sh data/intermediate_files/45-samplelist.txt data/intermediate_files/02-bis_set_real_unrelset data/intermediate_files/45-pcair_43noneur_unrelset

samplefile=$1
plinkfile=$2
outfile=$3

#tmpnum=21080301
tmpnum=$(date '+%Y%m%d%H%M%S')
mkdir "$tmpnum"
tmpfile1="${tmpnum}/1"
echo "Subsetting samples"
plink2 --bfile "$plinkfile" --keep "$samplefile" --make-bed --out "$tmpfile1" --memory 10000

infile="$tmpfile1"
tmpfile2="${tmpnum}/2"
echo "Preparing file for PCA"
bash ~/scripts/prePCA_plink2.sh "$infile" "$tmpfile2"

infile="${tmpfile2}-pruned.bed"
prefix="${tmpnum}/3-kinship"
tmpfile3="${tmpnum}/3.fam"
echo "Running KING to obtain kinship matrix"
awk '{print $1,$2,$3,$4,$5,"-9"}' "${infile%.bed}.fam" > $tmpfile3
king -b $infile --kinship --prefix $prefix --fam $tmpfile3

infile="${infile%.bed}"
kinmat="$prefix"
Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$infile" "$kinmat" "$outfile"

echo ""
echo "Checking what is the highest relationship pair"
[ -f "${prefix}.kin0" ] && maxkinship=$(sed '1d' ${prefix}.kin0 |awk -v max=0 '{if($8>max){want=$8; max=$8}}END{print want} ')
echo "Highest kinship found: ${maxkinship}"
awk -v maxkinship=$maxkinship '{ 
  if (maxkinship > 0.354) print "Duplicate/MZ twin found";
  else if (maxkinship >= 0.177 && maxkinship <= 0.354) print "1st-degree found";
  else if (maxkinship >= 0.0884 && maxkinship <= 0.177) print "2nd-degree found";
  else if (maxkinship >= 0.0442 && maxkinship <= 0.0884) print "3rd-degree found";
  else print "Unknown relationship of ",maxkinship," found";
}' <(sed '1d' ${prefix}.kin0 |awk -v max=0 '{if($8>max){want=$8; max=$8}}END{print want} ')


rm -r "${tmpnum}"

