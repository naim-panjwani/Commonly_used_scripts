# Script to convert PLINK binary format to VCF format
# Syntax: bash plink_to_vcf.sh <binary PLINK filename> <output filename without vcf extension>
# Example: bash plink_to_vcf.sh 04_omni_QCd_dup_rm_DI_rm 05_omni_QCd_dup_rm

echo "1. Generating frequency file using PLINK"
plink --noweb --bfile $1 --freq --out freq_temp >05_linkage_to_vcf.log

echo "2. Setting reference alleles"
awk '{print $2,"\t",$4}' freq_temp.frq |sed '1d' >05_ref_alleles_temp.txt
plink --noweb --bfile $1 --reference-allele ref_alleles_temp.txt --make-bed --out 05_ref_alleles_set_temp

echo "3. Initiating PLINK-SEQ new project"
pseq 05_temp_project new-project
echo "4. Loading PLINK file into PLINK-SEQ"
pseq 05_temp_project load-plink --file 05_ref_alleles_set_temp --id 05_ref_alleles_set_temp
echo "5. Writing VCF file using PLINK-SEQ"
pseq 05_temp_project write-vcf |gzip >$2.vcf.gz

