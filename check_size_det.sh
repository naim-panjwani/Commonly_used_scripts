#!/bin/bash
if [ ! -e check_size_det.out ]; then touch check_size_det.out; fi
cd fanwang/
for i in $(find -maxdepth 1 -type d -name '[!.]*'); do du -sh $i >>../check_size_det.out; done

# qsub check_size_det.sh -l walltime=167:59:00 -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -o ./jobout -e ./jobout -d `pwd` -N calcFanFolderSizeDet
