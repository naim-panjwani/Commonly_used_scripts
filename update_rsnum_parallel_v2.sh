#!/bin/bash
# Takes a bim file and a bed-formatted tabix indexed file from dbSNP and outputs a file with
# the subset of rs snp id's corresponding to the base-pair positions in the bim file
# This version breaks down the job into <num_cores> threads and then re-combines
# Syntax: bash <bim filename> <dbsnp filename> <output filename> <num_cores>
# Example: bash update_rsnum_parallel_v2.sh 38_ceu_cts_IIDupdated_18nov2014update.bim snp141.bed.gz  38_new_bimfile 30

bim=$1
num_lines=`wc -l ${bim} |awk '{print $1}'`
num_cores=$4
line_split=`echo $((num_lines / num_cores))`
outname=$3
split -d -a 3 -l ${line_split} $1 ${bim}.tmp.list
i=0
while [ $i -le $num_cores ]; do
  nohup bash update_rsnum_v2.sh ${bim}.tmp.list`printf "%03d" $i` $2 ${outname}`printf "%03d" $i`.txt >${outname}`printf "%03d" $i`.out &
  let i=$i+1
done
echo "Jobs sent. When done, execute combine_and_cleanup.sh script"
#touch ${outname}
#cat ${outname}*txt >${outname}
