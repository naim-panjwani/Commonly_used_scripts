#!/bin/bash
if [ ! -e check_size.out ]; then touch check_size.out; fi
cd fanwang/
find -maxdepth 1 -type d -name '[!.]*' -exec du -sh {} + >>../check_size.out

# qsub check_size.sh -l walltime=167:59:00 -l nodes=1:ppn=1 -l mem=8g -l vmem=8g -o ./jobout -e ./jobout -d `pwd` -N calcFanFolderSize
