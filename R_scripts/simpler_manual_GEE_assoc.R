setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/150615-UKRE_associations/CTS/")

plink_filename <- "05_CTS_SS_elp4_wout_Sardinians"
rawfilename <- "06_CTS_SS_elp4_wout_Sardinians.raw"
pcafilename <- "kingpcapc_no_sardinians.ped"
pca_rm_list <- "pca_exclude.txt"        # if no file, just enter the empty string ""
outname <- "07_assoc_gee_indep_wout_Sardinians.txt"

library(geepack)

pca <- read.table(pcafilename)
if(length(which(is.na(pca[,7])))>0) pca <- pca[-which(is.na(pca[,7])),]  # remove all NA lines
if(!(pca_rm_list %in% c("NA","0", ""))) {  # if a file is provided for exclusion of individuals in the PCA
  pca_remove_list <- read.table(pca_rm_list)
  pca <- subset(pca, !(as.character(pca[,2]) %in% as.character(pca_remove_list[,2]))) # individuals listed for exclusion are removed from the PCA file
}

fam <- read.table(paste0(plink_filename,".fam"), header=F)
fam <- cbind(index=1:nrow(fam), fam)
newfam <- fam[order(fam[,2]),]  # ordering is important for GEE to work properly as clusters (family members) need to be in adjacent rows
newpca <- pca # just to stick with the filename convention


for(i in 1:nrow(newfam)) { # the PCA file individuals must be in the same order as newfam
  nextid <- as.character(newfam[,3][i])
  newpca[i,] <- subset(pca, as.character(pca[,2]) %in% nextid)
}

newpca2 <- newpca
fam2 <- newfam[,-1]

# newpca2 and fam2 have the same order


# 1. Load data
map <- read.table(paste0(plink_filename,".bim"), header=F)
num.individuals<-nrow(fam2)
num.snps <- nrow(map)

print("Loading raw file")
ped <- read.table(rawfilename, header=T)  # May take awhile... not recommended for huge datasets; use bigmemory implementation for big datasets


# 2. Order the data
print("Ordering data")
ped <- ped[as.numeric(newfam[,1]),]
newped2 <- as.matrix(ped[,7:ncol(ped)])  # matrix takes less memory than data.frame and should speed things up
print("Data re-ordering done")


# Order of fam2, newpca2, ped and newped2 should all match at this point


# 3. Generate covariates variable
print("Generating covariates")
covar <- NULL
covar <- cbind(fam2[,1:2]) # FID/IID
sex <- fam2[,5]
PCs <- newpca2[,7:9] # 3 PCs -- MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
covar <- cbind(covar, sex, PCs)
names(covar) <- c("FID", "IID", "SEX", "PC1", "PC2","PC3")  # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED

snpnames <- as.character(map$V2)
chr <- map[,1]
CTS <- as.numeric(fam2[,6]) - 1
FIDs <- as.numeric(fam2[,1])  # for this to work, class(fam2[,1]) %in% "factor" must be TRUE


# Order of covar should match with newped2 and all others (fam2, newpca2, and ped)


print("Running GEE independence association model")
assoc.wrapper <- function(Genotype, snp.chr, snp.pos, snp.name, pheno, covariates, clusters) {
  SEX <- as.numeric(covariates$SEX) -1
  PC1 <- covariates$PC1    # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
  PC2 <- covariates$PC2
  PC3 <- covariates$PC3
  datum <- data.frame(cts=pheno, SEX=SEX, Genotype=Genotype, PC1=PC1, PC2=PC2, PC3=PC3, clusters=clusters)      # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
  datum <- na.omit(datum) # must omit missing data as GEE will fail otherwise
  tryCatch( # to deal with the situation where a particular GEE model does not converge or gives an error; the return will be NA
{ # first part of tryCatch -- executes first
  model <- geeglm(cts ~ SEX + Genotype + PC1 + PC2 + PC3, family="binomial",    # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
                  data=datum, corstr="independence", id=clusters)
  coef.table <- coef(summary(model))
  wald <- as.numeric(coef.table$Wald[which(rownames(coef.table) == "Genotype")])
  return(c(CHR=as.character(snp.chr), # character type so that the entire vector is character type
           POS=as.character(snp.pos),
           SNP=as.character(snp.name),
           beta=as.numeric(coef.table$Estimate[which(rownames(coef.table) == "Genotype")]),
           SE=as.numeric(coef.table$Std.err[which(rownames(coef.table) == "Genotype")]),
           OR=exp(as.numeric(coef.table$Estimate[which(rownames(coef.table) == "Genotype")])),
           Wald=wald,
           P.value=pchisq(wald,df=1,lower.tail=F))) 
} # if an error occurs, then the line below executes and that is what is returned instead
, error = function(e) {return(c(CHR=as.character(snp.chr), POS=as.character(snp.pos), SNP=as.character(snp.name), beta=NA, SE=NA, OR=NA, Wald=NA, P.value=NA))}
  ) 
}



SNPs_of_interest <- c("rs1495855", "rs17014744", "rs662702")

snp_indexes <- which(as.character(map[,2]) %in% SNPs_of_interest) # if the SNP is not found, it won't tell you and will not report!

results <- NULL
for(i in snp_indexes) {
  association <- assoc.wrapper(Genotype=as.numeric(newped2[,i]),
                               snp.chr=as.numeric(map[i,1]),
                               snp.pos=as.character(map$V4[i]),
                               snp.name=as.character(map$V2[i]),
                               pheno=CTS,
                               covariates=covar,
                               clusters=FIDs)
  results <- rbind(results, association)
}


print("Writing results to file")
write.table(results,outname, quote=F, row.names=F, col.names=T)



# Also do regular logistic regression - saving as a big list this time
lmmodels <- list()
length(lmmodels) <- length(snp_indexes)
names(lmmodels) <- map[snp_indexes,2]
counter<-0
for(i in snp_indexes) {
  counter <- counter+1
  datum <- data.frame(cts=CTS, SEX=covar$SEX, Genotype=as.numeric(newped2[,i]), PC1=covar$PC1, PC2=covar$PC2, PC3=covar$PC3, clusters=FIDs)  # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
  lmmodels[[counter]] <- glm(cts ~ SEX + Genotype + PC1 + PC2 + PC3, family="binomial", data=datum)   # MODIFY ACCORDING TO NUMBER OF PC's REQUIRED
}
lapply(lapply(lmmodels, summary), coef)

