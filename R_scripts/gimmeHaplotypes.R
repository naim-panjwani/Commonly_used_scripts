#!/usr/bin/R
# Input: Compressed VCF file with phased genotypes and snps_of_interest.txt file (one column of desired SNPs to analyze)
# Output: Four tab-delimited files are generated prefixed with <out_filename>: 
#         1) _hapMapfile.txt
#         2) _haplotypes_per_individual.txt
#         3) _haplotype_stats.txt
#         4) _diHaplotypes_stats.txt
# Syntax: Rscript gimmeHaplotypes.R <compressed_vcf_file> <snps_to_assess.txt> <out_filename_prefix>
# Example: Rscript gimmeHaplotypes.R 08_UKomni25for702_dups_nameUpdate_chr11.imputed.snps_of_interest.vcf.gz \
#                                    snps_of_interest.txt \
#                                    10_UKomni25for702


args=(commandArgs(TRUE))
vcf_filename <- as.character(args[1])
snps_to_assess_filename <- as.character(args[2])
out_prefix <- as.character(args[3])

library(data.table)
vcf <- fread(paste("gunzip -c ", vcf_filename), stringsAsFactors = F)
snps_to_assess <- fread(snps_to_assess_filename, stringsAsFactors=F, header=F)$V1
vcf <- subset(vcf, vcf$ID %in% snps_to_assess)

# Ensure first field is GT field and data is phased:
stopifnot(substring(vcf$FORMAT[1], 1,2) == "GT" & substring(vcf[1,10], 2,2) == "|")

all_haplotypes1 <- NULL
all_haplotypes2 <- NULL
for(i in 10:ncol(vcf)) {  # for each individual
  # Get the haplotype among the snps to assess
  haplotype1 <- ""
  haplotype2 <- ""
  for(j in 1:nrow(vcf)) haplotype1 <- paste0(haplotype1, as.numeric(substring(vcf[j,i,with=F],1,1)))
  for(j in 1:nrow(vcf)) haplotype2 <- paste0(haplotype2, as.numeric(substring(vcf[j,i,with=F],3,3)))
  
  #Store it as a string:
  all_haplotypes1 <- c(all_haplotypes1, haplotype1)
  all_haplotypes2 <- c(all_haplotypes2, haplotype2)
}
all_haplotypes <- c(all_haplotypes1, all_haplotypes2)
mapfile <- vcf[,c(1:5,8),with=F]
haps_per_ind <- data.frame(IID=colnames(vcf)[10:ncol(vcf)], Haplotype1=all_haplotypes1, Haplotype2=all_haplotypes2)

x <- as.data.frame(cbind(summary(as.factor(all_haplotypes)), summary(as.factor(all_haplotypes))/length(all_haplotypes)*100))
x <- data.frame(Haplotype=rownames(x), N=x[,1], Percentage=x[,2])
x <- x[order(x[,3], decreasing = T),]
x <- cbind(x, HapCode=paste0("H",seq(1:nrow(x))))
hap_stats <- x

haps_per_ind2 <- NULL
for(i in 1:nrow(haps_per_ind)) {
  major_hap1 <- NA
  major_hap2 <- NA
  major_hap1 <- ifelse(length(which(as.character(x[,1]) %in% haps_per_ind[i,2]))==0, "Other", 
                 paste0("H", which(as.character(x[,1]) %in% haps_per_ind[i,2])))
  major_hap2 <- ifelse(length(which(as.character(x[,1]) %in% haps_per_ind[i,3]))==0, "Other", 
                       paste0("H", which(as.character(x[,1]) %in% haps_per_ind[i,3])))
  haps_per_ind2 <- rbind(haps_per_ind2, c(IID=as.character(haps_per_ind[i, 1]), 
                                                      Haplotype1=as.character(haps_per_ind[i, 2]), 
                                                      Haplotype2=as.character(haps_per_ind[i, 3]),  
                                                      Major1=major_hap1, Major2=major_hap2))
}

diHaps <- NULL
for(i in 1:nrow(haps_per_ind2)) diHaps <- c(diHaps, paste0(haps_per_ind2[i,4], haps_per_ind2[i,5]))
diHaplos <- as.data.frame(summary(as.factor(diHaps))/length(diHaps)*100)
diHaplos <- data.frame(DiHapCode=rownames(diHaplos), N=summary(as.factor(diHaps)), Percentage=summary(as.factor(diHaps))/length(diHaps)*100)
diHaplos <- diHaplos[order(diHaplos[,3], decreasing=T),]

haps_per_ind3 <- cbind(haps_per_ind2, diHaps)


write.table(mapfile, paste0(out_prefix, "_mapfile.txt"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(haps_per_ind3, paste0(out_prefix, "_haplotypes_per_individual.txt"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(hap_stats, paste0(out_prefix, "_haplotype_stats.txt"), quote=F, row.names=F, col.names=T, sep="\t")
write.table(diHaplos, paste0(out_prefix, "_diHaplotypes_stats.txt"), quote=F, row.names=F, col.names=T, sep="\t")


