#!/bin/bash
# Split a PLINK file into #cores chunks and send #cores jobs and then recombine results
# Dependency: calculate_num_splits.R
# Syntax: bash parallel_kinship.sh <PLINK_filename> <output_filename> <num_cores>
# Example: bash parallel_kinship.sh 09_pruned 15_kinship 30 

PLINK_filename=$1
output_filename=$2
num_cores=$3
freq_file=${PLINK_filename}_freq
num_ind=`wc -l ${PLINK_filename}.fam |cut -d' ' -f1`

# Generate frequency file:
plink --noweb --bfile ${PLINK_filename} --freq --out ${freq_file}

# Split files into $num_cores chunks
Rscript calculate_num_splits.R ${num_cores}
num_splits=`cat num_splits.txt`
num_ind_per_file=`echo $((num_ind / num_splits))`
#let num_splits=$num_splits-1
gawk '{print $1,$2}' ${PLINK_filename}.fam | split -d -a 3 -l $num_ind_per_file - tmp.list

# Send the pruning jobs
let i=0 a=${num_splits}
  let j=0
    while [ $i -le $a ]
      do
        while [ $j -le $a ]
        do
          nohup plink --noweb --bfile $PLINK_filename \
                 --read-freq ${freq_file}.frq \
                 --genome \
                 --genome-lists tmp.list`printf "%03i\n" $i` \
                                tmp.list`printf "%03i\n" $j` \
                 --min 0.05 \
                 --out ${output_filename}.sub.$i.$j &
        let j=$j+1
        done
      let i=$i+1
      let j=$i
    done
let i=0
let j=0
PATH=$PATH:/home/naim/scripts/
while [ $i -le $a ]
do
  while [ $j -le $a ]
  do
    waiton $(findpid sub.${i}.${j})
    let j=$j+1
  done
let i=$i+1
let j=$i
done

files=`ls ${output_filename}.sub*genome`
(head -1 ${output_filename}.sub.0.0.genome; for i in $files; do sed '1d' $i; done) >${output_filename}.genome

#rm tmp.list*
#rm ${output_filename}.sub.*
