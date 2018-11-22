#!/bin/bash

file=$1

while read line; do
  chr=$(echo $line |cut -d' ' -f1)
  startBP=$(echo $line |cut -d' ' -f2)
  endBP=$(echo $line |cut -d' ' -f3)
  (gunzip -c $file |head -1; gunzip -c $file |awk -v chr=$chr -v start=$startBP -v end=$endBP '$1==chr && $2>=start && $2<=end') |bgzip -c >subsets/${file%.txt.gz}_chr${chr}_${startBP}_${endBP}.txt.gz
done < regions_of_interest.txt

