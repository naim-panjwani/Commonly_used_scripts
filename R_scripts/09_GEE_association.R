# GEE
# CTS = Genotype + Sex + Site + PCs

library(geepack)

# 1. Load data (convert to ped/map first)
#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/141020-GWAS_GEE_TOR_CF_Ctrls_CTS_exome")
ped <- read.table("newpped.ped", header=F)
map <- read.table("08_chr22_merged_common_snps_only_recoded.map", header=F)
newped <- ped[order(ped$V1),] #re-order ped file so that same FID (cluster) is beside each other


# 2. Load PCs
#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/141020-GWAS_GEE_TOR_CF_Ctrls_CTS_exome/pca/")
pca <- read.table("./pca/15_mergedset_pca_outliers_rm_with_relatedspc.ped")
newpca <- pca[order(pca$V1),1:10]

# 3. Remove the 7 PCA outliers (all were controls)
newped2 <- subset(newped, as.character(newped$V2) %in% as.character(newpca$V2))

# 4. Load Site information
#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/")
RE_pheno_data <- read.csv("RE_merged_data_correction1.csv")
RE_pheno_data <- subset(RE_pheno_data, as.character(RE_pheno_data$IID) %in% as.character(newped2$V2))

table(RE_pheno_data$Site)
#Argentina    Canada    France    Kerala  Sardinia        UK       USA 
#0            13        6         0       0               46       62 
# 1061_301 got left behind because forgot to do the IID corrections on imputed data (formerly tagged as 1061_302)


# 5. Add the control individuals to RE_pheno_data with dummy variables, and put both RE_pheno_data and newpca in equal orders
temp <- subset(newped2[,2], newped2$V6==1)
newtemp <- matrix(nrow=length(temp),ncol=dim(RE_pheno_data)[2])
newtemp[,1] <- as.character(temp)
colnames(newtemp) <- colnames(RE_pheno_data)
#new_RE_pheno_data <- cbind(IID=as.character(RE_pheno_data$IID), RE_pheno_data[,2:dim(RE_pheno_data)[2]])
new_RE_pheno_data <- rbind(RE_pheno_data,newtemp)
# Put it in order:
new_RE_pheno_data2 <- new_RE_pheno_data
newpca2 <- subset(newpca, as.character(newpca$V2) %in% as.character(newped2$V2))
for(i in 1:dim(new_RE_pheno_data)[1]) {
  nextid <- as.character(newped2$V2[i])
  new_RE_pheno_data2[i,] <- subset(new_RE_pheno_data, as.character(new_RE_pheno_data$IID) %in% nextid)
  newpca2[i,] <- subset(newpca, as.character(newpca$V2) %in% nextid)
}


# 6. Ensure order of IID's is the same in newped2 and RE_pheno_data
is_in_order=TRUE
not_in_order_indices=NULL
print("Is new_RE_pheno_data2 in order?")
for(i in 1:dim(newped2)[1]) {
  if(as.character(newped2$V2[i]) != as.character(new_RE_pheno_data2$IID[i])) {
    is_in_order=FALSE
    not_in_order_indices=c(not_in_order_indices,i)
    print(paste(i, newped2$V2[i], new_RE_pheno_data2$IID[i]))
  }
}
print(paste("is_in_order =",is_in_order))

is_in_order=TRUE
not_in_order_indices=NULL
print("Is newpca2 in order?")
for(i in 1:dim(newped2)[1]) {
  if(as.character(newped2$V2[i]) != as.character(newpca2$V2[i])) {
    is_in_order=FALSE
    not_in_order_indices=c(not_in_order_indices,i)
    print(paste(i, newped2$V2[i], newpca2$V2[i]))
  }
}
print(paste("is_in_order =",is_in_order))


# 7. Generate covariates variable and file
covar <- cbind(newped2[,1:2]) # FID/IID
sex <- newped2$V5
site <- as.character(new_RE_pheno_data2$Site)
site[is.na(site)] <- "TOR_CF"
site <- as.factor(site)
PC4 <- newpca2[,7:10]

covar <- cbind(covar, sex, site, PC4)
names(covar) <- c("FID", "IID", "sex", "site", "PC1", "PC2", "PC3", "PC4")
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/141020-GWAS_GEE_TOR_CF_Ctrls_CTS_exome")
write.table(covar, "09_covariates.txt",quote=F,row.names=F,col.names=T)
write.csv(covar, "09_covariates.csv",quote=F,row.names=F)


# 8. Model the genotypes for each SNP
source("AssocCG.R")
# MinorAllele(newped2[,7],newped2[,8])
# GenCounts(newped2[,7],newped2[,8])
# as.numeric(Genotype(newped2[,7],newped2[,8], AA_code=2, Aa_code=1, aa_code=0)$Genotype.code)
# SNPDatum(newped2, newped2[,7],newped2[,8])
# analysis <- Table23(newped2[,1:10], map)

Convenient.Genotype <- function(allele1, allele2) {
  return(as.numeric(Genotype(allele1, allele2, AA_code=2, Aa_code=1, aa_code=0)$Genotype.code))
}

num.individuals <- dim(newped2)[1]
num.snps <- dim(map)[1]
geno_data <- matrix(nrow=num.individuals, ncol=num.snps)
snpnames <- as.character(map$V2)

library(doMC)
library(foreach)
registerDoMC(60)

# ptm <- proc.time()
# for(snp in 1:15) {
#   col.position <- 7 + 2*(snp-1)
#   geno_data[,snp] <- Convenient.Genotype(newped2[,col.position], newped2[,col.position+1])
# }
# proc.time() - ptm

ptm <- proc.time()
geno_data <- foreach(snp=1:num.snps, .combine=cbind) %dopar% {
  col.position <- 7 + 2*(snp-1)
  Convenient.Genotype(newped2[,col.position], newped2[,col.position+1])
}
proc.time() - ptm
# colnames(geno_data) <- snpnames
rownames(geno_data) <- as.character(newped2[,2])
# write.table(geno_data,"genotype_data.txt",quote=F,row.names=T,col.names=T)
# write.csv(geno_data,"genotype_data.txt",quote=F,row.names=T)
# write.csv(geno_data,"genotype_data.txt",quote=F,row.names=F)

# 9. Model the association (GEE independence correlation structure)
model <- list()
length(model) <- num.snps
cts <- newped2[,6] - 1
SEX <- as.factor(covar$sex - 1)
SITE <- as.factor(covar$site)
for(i in 1:15) {
  snpname <- as.character(map$V2[i])
  model[[i]] <- geeglm(cts ~ SEX + geno_data[,i] + SITE + covar$PC1 + covar$PC2 + covar$PC3 + covar$PC4,
                         family="binomial", corstr="independence", id=newped2[,1])
  names(model)[i] <- snpname
}
summary(model[[1]])

