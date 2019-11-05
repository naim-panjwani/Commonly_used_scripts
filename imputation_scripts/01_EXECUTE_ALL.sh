echo "Renaming files"
bash 01_rename_files.sh

echo "1. Manually creating list of twins/duplicates to remove from dataset"

echo "4007_303 4007_303
1030_302 1030_302
LKS_I LKS_I
S165 S165
7042_302 7042_302
7072_301 7072_301" >02_twins_duplicates_to_remove.txt

echo "2. Removing twins/duplicates with PLINK"
plink --noweb --bfile 10_nof_rm --remove 02_twins_duplicates_to_remove.txt --make-bed --out 02_omni_QCd_dup_rm >02_remove_duplicates.log

echo "3. Running ped_build_REversion3.R to generate updated fam file to include trio relationships"
R CMD BATCH 03_ped_build_REversion3.R

echo "4. Updating 02_omni_QCd_dup_rm.fam file"
cp 03_pedBuildOmniChip.fam 02_omni_QCd_dup_rm.fam

echo "5. Removing indels"
awk '!($5 == "A" || $5 == "T" || $5 == "G" || $5 == "C") || !($6 == "A" || $6 == "T" || $6 == "G" || $6 == "C")' 02_omni_QCd_dup_rm.bim |cut -f2 >DIsnps.txt
plink --noweb --bfile 02_omni_QCd_dup_rm --exclude DIsnps.txt --make-bed --out 04_omni_QCd_dup_rm_DI_rm

echo "5. Converting PLINK format files to VCF format using PLINK-SEQ"
bash 05_plink_to_vcf.sh 02_omni_QCd_dup_rm 05_omni_QCd_dup_rm_DI_rm

echo "6. Generating non_PASSing_snps.txt containing the list of SNPs that did not get PASS in VCF file"
gunzip -c 05_omni_QCd_dup_rm.vcf.gz |grep -v ^# |grep -v PASS |cut -f3 >non_PASSing_snps.txt

echo "7. Renaming the file for convenience:"
cp 05_omni_QCd_dup_rm.vcf.gz 06_Omni.vcf.gz

echo "8. Preparing sample file for conform-gt.jar program"
echo "8.1 Removing chr prefix"
gunzip -c 06_Omni.vcf.gz |grep -v 'rs10005853' |sed 's/^0/X/g' |bgzip >07_Omni_chrXrepaired.vcf.gz
