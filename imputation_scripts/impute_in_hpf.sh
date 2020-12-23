#!/bin/bash

for chri in 1; do
chr=$chri
if [[ $chri == "23" ]]; then chr="X"; fi
blockcount=`ls 13_chr${chr}_UKRE_Spit_for_Sci.*.vcf.gz |wc -l`
for i in `ls 13_chr${chr}_UKRE_Spit_for_Sci.*.vcf.gz`
do
block=`echo ${i} |cut -d'.' -f2`
if [[ $block == 1 ]]; then 
  minbp=`zcat 1000Genomes_phase3_for_BEAGLE_v4/chr${chr}.1kg.phase3.v5.vcf.gz |grep -v ^'#' |head -1 |cut -f2`
  maxbp=`zcat 13_chr${chr}_UKRE_Spit_for_Sci.${block}.vcf.gz |grep -v ^'#' |tail -1 |cut -f2`
elif [[ $block == $blockcount ]]; then
  minbp=`zcat 13_chr${chr}_UKRE_Spit_for_Sci.${block}.vcf.gz |grep -v ^'#' |head -1 |cut -f2`
  maxbp=`zcat 1000Genomes_phase3_for_BEAGLE_v4/chr${chr}.1kg.phase3.v5.vcf.gz |grep -v ^'#' |tail -1 |cut -f2`
else
  minbp=`zcat 13_chr${chr}_UKRE_Spit_for_Sci.${block}.vcf.gz |grep -v ^'#' |head -1 |cut -f2`
  maxbp=`zcat 13_chr${chr}_UKRE_Spit_for_Sci.${block}.vcf.gz |grep -v ^'#' |tail -1 |cut -f2`
fi

qsub -l mem=120g,vmem=120g,nodes=1:ppn=32,walltime=720:00:00 \
     -o /home/naim/UKRE_imputation/150203_UKRE_plus_ALL_Spit_for_Science/joboutdir/ \
     -e /home/naim/UKRE_imputation/150203_UKRE_plus_ALL_Spit_for_Science/joboutdir/ \
     -d /home/naim/UKRE_imputation/150203_UKRE_plus_ALL_Spit_for_Science/ \
     -N chr${chr}.block${block} \
     -v chr=${chr},block=${block},minbp=${minbp},maxbp=${maxbp} \
       runjava.pbs
done
done

echo "End of script"

