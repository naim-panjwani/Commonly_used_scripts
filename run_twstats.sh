#!/bin/bash

module load gsl

eigenvalues=$1
outname=$2

if [ ! -e twtable ]; then
  ln -s ~/scripts/twtable
fi

numcols=$(head -1 $eigenvalues |wc -w)
if [[ $numcols -eq 1 ]]; then
  twstats -t twtable -i "$eigenvalues" -o "$outname"
elif [[ $numcols -eq 2 ]]; then
  echo "Assuming eigenvalues are listed on the 2nd column"
  cut -f2 "$eigenvalues" > tmp211216
  twstats -t twtable -i tmp211216 -o "$outname"
  rm tmp211216
else
  echo "Too many columns detected. Check your file"
  exit
fi

