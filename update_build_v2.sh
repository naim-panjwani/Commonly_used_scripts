#!/bin/sh

# This script has been fixed for the chr25 (PAR region of X). 
# The original script assummed chr25 labeled as chr23 instead.
# Also, the final step of removing SNPs not in the strand file is avoided
# as this seems to take out quite a few SNPs (about 1-3%)

#A script for updating a binary ped file using one of Will's strand files
#NRR 17th Jan 2012

#V2 13th Feb 2012. Added code to retain only SNPs in the strand file

#Required parameters:
#1. The original bed stem (not including file extension suffix)
#2. The strand file to apply
#3. The new stem for output
#Result: A new bed file (etc) using the new stem

#Unpack the parameters into labelled variables
stem=$1
strand_file=$2
outstem=$3
echo Input stem is $stem
echo Strand file is $strand_file
echo Output stem is $outstem


#Cut the strand file into a series of Plink slices
chr_file=$strand_file.chr
pos_file=$strand_file.pos
flip_file=$strand_file.flip
cat $strand_file | cut -f 1,2 > $chr_file
cat $strand_file | cut -f 1,3 > $pos_file
cat $strand_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file

#Because Plink only allows you to update one attribute at a time, we need lots of temp
#Plink files
temp_prefix=TEMP_FILE_XX72262628_
temp1=$temp_prefix"1"
temp2=$temp_prefix"2"
temp3=$temp_prefix"3"
temp4=$temp_prefix"4"
temp5=$temp_prefix"5"
temp6=$temp_prefix"6"

# Naim custom pre-steps:
# Need to update chr "25" (the PAR region of X) to 23
grep ^"25\s" ${stem}.bim |cut -f2,1 |sed 's/^25/23/g' |awk '{print $2"\t"$1}' > $temp1
plink --allow-no-sex --bfile $stem --update-chr $temp1 --make-bed --out $temp2


#1. Apply the chr
plink --allow-no-sex --bfile $temp2 --update-map $chr_file --update-chr --make-bed --out $temp3
#2. Apply the pos
plink --allow-no-sex --bfile $temp3 --update-map $pos_file --make-bed --out $temp4
#3. Apply the flip
plink --allow-no-sex --bfile $temp4 --flip $flip_file --make-bed --out $temp5
#4. Re-assign chr25
awk '{print $2"\t"$1}' $temp1 |sed 's/^23/25/g' |awk '{print $2"\t"$1}' > $temp6
plink --allow-no-sex --bfile $temp5 --update-chr $temp6 --make-bed --out $outstem


#4. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file <-- Naim: disabling this step as it seems to rid of many SNPs...
#plink --allow-no-sex --bfile $temp5 --extract $pos_file --make-bed --out $outstem # <-- disabling this step

#Now delete any temporary artefacts produced
rm -f $temp_prefix*

