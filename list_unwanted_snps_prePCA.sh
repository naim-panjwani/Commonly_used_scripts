#!/bin/bash
# Syntax: bash list_unwanted_snps_prePCA.sh <bim filename> <output filename>

grep ^23 $1 |cut -f2 >$2
awk '$5 == "I" || $6 == "I" || $5 == "D" || $6 == "D"' $1 |cut -f2 >>$2
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' $1 |cut -f2 >>$2

