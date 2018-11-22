#!/bin/bash

dataset=$1 #ensure to run prePCA.sh on it first
prefix=$2


echo "Starting KING PCA"
king -b ${dataset}.bed --pca 20 --prefix ${prefix}_mergedset_pca >${prefix}_king_pca.out
echo "KING PCA DONE"
echo "Generating Eigenvalues file"
grep 'eigenvalues' ${prefix}_king_pca.out |cut -d' ' -f4- |tr ' ' '\n' >${prefix}_eigenvalues.txt
echo "Tracy Widom Analysis"
bash ~/scripts/tracy-widom.sh ${prefix}_eigenvalues.txt ${prefix}_tracy-widom.txt
echo "Replacing X's in PCA results file by NA's"
sed 's/X/NA/g' ${prefix}_mergedset_pcapc.ped >${prefix}tmp
cat ${prefix}tmp >${prefix}_mergedset_pcapc.ped
rm ${prefix}tmp

