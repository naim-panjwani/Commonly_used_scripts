#!/bin/bash

infile=$1

tmpnum="210818"

module load vep/94
dbsnp="data/raw_data/dbSNP151_GRCh38p7.vcf.gz"
annotated="${infile%.txt}_annotated.txt"
cachedir="/hpf/largeprojects/struglis/datasets/vep_cache"
tmpfile="tmp${tmpnum}"
outfile="${infile%.txt}_rsids.txt"
sed 's/^23/X/g' $infile > $tmpfile
#cp $tmpfile $infile
vep -i $tmpfile -o $annotated --cache --offline --dir_cache $cachedir --check_existing --nearest symbol --fork 2 --force_overwrite

(head -1 "$tmpfile"; LANG=en_EN join -1 3 -2 1 -o 1.1 1.2 2.2 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 <(LANG=en_EN sort -k3 "$tmpfile") <(grep -v ^"##" "$annotated" |cut -f1,13 |LANG=en_EN sort -k1) |uniq |tr ' ' '\t' |LANG=en_EN sort -k2 -n) > "$outfile"

