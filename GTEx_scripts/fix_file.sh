#!/bin/bash

file=$1

paste <(echo -e "chr\tpos"; gunzip -c $file |sed '1d' |cut -f2 |tr '_' '\t' |awk '{print $1"\t"$2}') <(gunzip -c $file) |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$6"\t"$7"\t"$8}' |bgzip -c >${file%.txt.gz}_fixed.txt.gz


