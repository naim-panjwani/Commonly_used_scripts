#!/bin/bash

INFO_column=$1
# writes to stdout
# ensure HEADER is not present
# INFO column format MUST BE in the following order: AR2=0.00;DR2=0.00;AF=0.0020;IMP (ie. AR2 and AF must be in the 1st and 3rd positions of the semicolon-delimited INFO field)

cat ${INFO_column} |tr ';' '\t' |tr '=' '\t' |cut -f2,6 >temp171221
while read line; do if [[ $line == *","* ]]; then echo $line |cut -f2 |bash get_largest_multiallelic.sh -; else echo $line |cut -f2; fi <temp171221 
(echo -e "Rsq\tALT_Frq\tMAF"; paste temp171221 <(awk '{if($2>0.5) {print 1-$2} else {print $2}}' temp171221)) 
rm temp171221 


