#!/hpf/tools/centos6/R/3.1.1/bin/R

args=(commandArgs(TRUE))

bplink_filename <- as.character(args[1])
rawplink_filename <- as.character(args[2])
pca_filename <- as.character(args[3])
output_filename <- as.character(args[4])

library(iterators)
library(doMC)
library(foreach)
registerDoMC(32)
library(bigmemory)

(num_cores_to_use <- getDoParWorkers())

#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/150320-19_Reinthaler_GABA_genes")
pca <- read.table(pca_filename)
if(length(which(is.na(pca[,7])))>0) pca <- pca[-which(is.na(pca[,7])),]

fam <- read.table(paste0(bplink_filename,".fam"), header=F)
fam <- cbind(index=1:nrow(fam), fam)
ctrls <- subset(fam, fam[,7]==1)
cases <- subset(fam, fam[,7]==2)
ctrls[,2] <- ctrls[,3]
newfam <- rbind(cases,ctrls)
newfam <- newfam[order(newfam[,2]),]
newpca <- pca
for(i in 1:nrow(newfam)) {
  nextid <- as.character(newfam[,3][i])
  newpca[i,] <- subset(pca, as.character(pca[,2]) %in% nextid)
}

newpca2 <- newpca
fam2 <- newfam[,-1]




# Load dataset
map <- read.table(paste0(bplink_filename,".bim"),header=F)
num.individuals<-dim(fam2)[1]
num.snps <- dim(map)[1]
print("Loading raw file")

ped <- read.big.matrix(paste0(rawplink_filename,".raw"),
                       type="short",
                       header=T,
                       sep=" ")
# 1. Order the data
print(paste("Ordering data - important if running GEE"))
newped <- cbind(fam[,-1], ped[,7:ncol(ped)])
newped2 <- newped[as.integer(newfam[,1]),] #re-order ped file so that same FID (cluster) is beside each other (gives problems with GEE otherwise)

options(bigmemory.typecast.warning=FALSE)
updated_ped <- as.big.matrix(as.matrix(newped2[,7:ncol(newped2)])
                             ,type="short"
)
newped2 <- updated_ped
print("Data re-ordering done")


# 3. Generate covariates variables
print(paste("Generating covariates"))
covar <- NULL
covar <- cbind(fam2[,1:2]) # FID/IID
sex <- fam2[,5]
PCs <- newpca2[,7:9] # 3 PCs 
covar <- cbind(covar, sex, PCs)
names(covar) <- c("FID", "IID", "SEX", "PC1", "PC2","PC3")
snpnames <- as.character(map$V2)

CTS <- as.numeric(fam2[,6]) - 1
FIDs <- as.numeric(factor(fam2[,1]))

# 4. Association
print("Base model summary without genotype: CTS|SSD = Sex + PC1 + PC2 + PC3")
Sex <- as.numeric(covar$SEX) -1
pc1 <- covar$PC1
pc2 <- covar$PC2
pc3 <- covar$PC3
base_model <- glm(CTS ~ Sex + pc1 + pc2 + pc3, family="binomial")
summary(base_model)


assoc.wrapper_full <- function(Genotype, snp.name, snp.pos, pheno, covariates, clusters) {
  SEX <- as.numeric(covariates$SEX) -1
  PC1 <- covariates$PC1
  PC2 <- covariates$PC2
  PC3 <- covariates$PC3
  cts <- pheno
  datum <- data.frame(cts=cts, SEX=SEX, Genotype=Genotype, PC1=PC1, PC2=PC2, PC3=PC3, clusters=clusters)
  datum <- na.omit(datum)
  model <- glm(cts ~ SEX + Genotype + PC1 + PC2 + PC3, family="binomial", data=datum)
  model_no_pcs <- glm(cts ~ SEX + Genotype, family="binomial", data=datum)
  coef.table <- coef(summary(model))
  zvalue <- as.numeric(coef.table[,3][which(rownames(coef.table) == "Genotype")])
  coef.table_noPC <- coef(summary(model_no_pcs))
  zvalue_noPC <- as.numeric(coef.table_noPC[,3][which(rownames(coef.table_noPC) == "Genotype")])
  diff_logLik=-2*(as.numeric(logLik(model_no_pcs)) - as.numeric(logLik(model)))
  return(c(SNP=as.character(snp.name),
           POS=as.character(snp.pos),
           beta_3PC=as.numeric(coef.table[,1][which(rownames(coef.table) == "Genotype")]),
           SE_3PC=as.numeric(coef.table[,2][which(rownames(coef.table) == "Genotype")]),
           OR_3PC=exp(as.numeric(coef.table[,1][which(rownames(coef.table) == "Genotype")])),
           Z_3PC=zvalue,
           P_3PC=2*pnorm(-abs(zvalue)),
           logLik_3PC=as.numeric(logLik(model)),
           beta_noPC=as.numeric(coef.table_noPC[,1][which(rownames(coef.table_noPC) == "Genotype")]),
           SE_noPC=as.numeric(coef.table_noPC[,2][which(rownames(coef.table_noPC) == "Genotype")]),
           OR_noPC=exp(as.numeric(coef.table_noPC[,1][which(rownames(coef.table_noPC) == "Genotype")])),
           Z_noPC=zvalue_noPC,
           P_noPC=2*pnorm(-abs(zvalue_noPC)),
           logLik_noPC=as.numeric(logLik(model_no_pcs)),
           diff_logLik=diff_logLik,
           diff_logLik_P=pchisq(diff_logLik,df=3,lower.tail=F))
         )
}


assoc.wrapper <- function(Genotype, pheno, covariates, clusters) {
  SEX <- as.numeric(covariates$SEX) -1
  PC1 <- covariates$PC1
  PC2 <- covariates$PC2
  PC3 <- covariates$PC3
  cts <- pheno
  datum <- data.frame(cts=cts, SEX=SEX, Genotype=Genotype, PC1=PC1, PC2=PC2, PC3=PC3, clusters=clusters)
  datum <- na.omit(datum)
  model <- glm(cts ~ SEX + Genotype + PC1 + PC2 + PC3, family="binomial", data=datum)
  coef.table <- coef(summary(model))
  zvalue <- as.numeric(coef.table[,3][which(rownames(coef.table) == "Genotype")])
  return(chi_3PC=zvalue^2)
}


print(paste("Running Simple Logistic association model"))
coefs <- foreach(i=1:num.snps, .combine=rbind) %dopar% {
  assoc.wrapper_full(Genotype=as.numeric(newped2[,i]), 
                     snp.name=as.character(map$V2[i]),
                     snp.pos=as.character(map$V4[i]),
                     pheno=CTS,
                     covariates=covar,
                     clusters=FIDs)  
}

# Get the observed sum stat
chi <- (as.numeric(coefs[,'Z_3PC']))^2
(chi_obs <- sum(chi))

# Save
print("Saving observed statistics")
write.table(coefs, output_filename, quote=F, row.names=F, col.names=T)

# Permutations
NN = 10000

# Shuffle phenotypes:
print(paste("Starting permutations N =",NN))
print("NOTE: TAKING GENOTYPE BIG.MATRIX INTO MEMORY AND PARALLELIZING BY PERMUTATION ITERATIONS")
print(paste0("This works best when the number of genotype markers (",num.snps, ") < number of permutations (",NN, ")"))

tmpped <- newped2[,]
null_chisquare_sums <- foreach(n=1:NN, .combine=c) %dopar% {
  sum(as.numeric(apply(tmpped,2,assoc.wrapper,pheno=sample(CTS),covariates=covar,clusters=FIDs)))
}

(empirical_p <- sum(null_chisquare_sums > chi_obs) / NN)
print(paste("Empirical P-value = ",empirical_p))


