setwd("/Users/naim panjwani/Documents/Strug/RE_full_impute")

vcf <- read.table("temp2")
sample <- read.table("sample_temp2")

s <- vcf$V6
allelic_r2 <- as.numeric(gsub("AR2=(.*?);.*","\\1",s))
allele_freq <- as.numeric(gsub(".*;AF=(.*?)","\\1",s))
summary(allelic_r2)
dim(vcf[which(allelic_r2>0.8),])
head(vcf[which(allelic_r2>0.8),])
hist(allelic_r2)

low_allelic_r2 <- allelic_r2[allelic_r2<0.2]
index_low_allelic_r2 <- which(allelic_r2<0.2)

plot(low_allelic_r2,allele_freq[index_low_allelic_r2])
hist(allele_freq[index_low_allelic_r2], main="Allele Frequency for\nSNPs with Allelic R2<0.2",
     xlab="Allele frequency")
# Most of these SNPs are rare variants


setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3/imputation/trial2/impute_chr22/")
setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3/imputation")

