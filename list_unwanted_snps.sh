#!/bin/bash
# Script to list all non-chromosome 1-23 SNPs and ambiguous A/T G/C SNPs
# Syntax: bash list_unwanted_snps.sh <bim filename> <output filename>

awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "G" && $6 == "C") || ($5 == "C" && $6 == "G")' $1 |cut -f2 >$2
grep -Ev '\b([1-9]|1[0-9]|2[0-3])\b' $1 >>$2

