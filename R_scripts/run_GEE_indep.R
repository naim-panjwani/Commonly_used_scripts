#!/hpf/tools/centos6/R/3.1.1/bin/R
# GEE Independence
# CTS = Genotype + Sex + PC1 + PC2 + PC3 + PC4
# Syntax: Rscript run_gee_indep.R <bPLINK_filename> <raw_filename> <pca_filename> <pca_exclude_list> <outfile_name>
# Example: Rscript run_GEE_indep.R UKRE_RDG_chr1 chr1.raw kingpcapc.ped king_pca_parents_with_missing_rdg.txt gee_exch_rdg_chr1_assoc.txt

args=(commandArgs(TRUE))

plink_filename <- as.character(args[1])
rawfilename <- as.character(args[2])
pcafilename <- as.character(args[3])
pca_rm_list <- as.character(args[4])
outname <- as.character(args[5])

library(doMC)
library(foreach)
registerDoMC(32)
library(bigmemory)
library(geepack)

(num_cores_to_use <- getDoParWorkers())


pca <- read.table(pcafilename)
if(length(which(is.na(pca[,7])))>0) pca <- pca[-which(is.na(pca[,7])),]
if(!(pca_rm_list %in% c("NA","0", ""))) {
  pca_remove_list <- read.table(pca_rm_list)
  pca <- subset(pca, !(as.character(pca[,2]) %in% as.character(pca_remove_list[,2])))
}

fam <- read.table(paste0(plink_filename,".fam"), header=F)
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




  # Load data
  map <- read.table(paste0(plink_filename,".bim"), header=F)
  num.individuals<-dim(fam2)[1]
  num.snps <- dim(map)[1]
  print("Loading raw file")
  ped <- read.big.matrix(rawfilename,
                         type="short",
                         header=T,
                         sep=" ")
  
  # 1. Order the data
  print("Ordering data")
#  newped <- cbind(fam[,-1], ped[,7:ncol(ped)])
#  newped2 <- newped[as.integer(newfam[,1]),] #re-order ped file so that same FID (cluster) is beside each other (gives problems with GEE otherwise)
  
#  options(bigmemory.typecast.warning=FALSE)
#  updated_ped <- as.big.matrix(as.matrix(newped2[,7:ncol(newped2)])
#                               ,type="short"
#                              )
mpermute(ped, order=as.numeric(newfam[,1]))
newped2 <- sub.big.matrix(ped, firstRow = 1, lastRow = nrow(ped), firstCol = 7, lastCol = ncol(ped))
   print("Data re-ordering done")
  
  
  # 3. Generate covariates variable
  print("Generating covariates")
  covar <- NULL
  covar <- cbind(fam2[,1:2]) # FID/IID
  sex <- fam2[,5]
  PCs <- newpca2[,7:10] # 4 PCs 
  covar <- cbind(covar, sex, PCs)
  names(covar) <- c("FID", "IID", "SEX", "PC1", "PC2","PC3", "PC4")
  snpnames <- as.character(map$V2)
  
  
  # 5. Model the association
  
  print("Running GEE independence association model")
  assoc.wrapper <- function(Genotype, snp.name, snp.pos, pheno, covariates, clusters) {
    SEX <- as.numeric(covar$SEX) -1
    PC1 <- covar$PC1
    PC2 <- covar$PC2
    PC3 <- covar$PC3
    PC4 <- covar$PC4
    cts <- pheno
    datum <- data.frame(cts=cts, SEX=SEX, Genotype=Genotype, PC1=PC1, PC2=PC2, PC3=PC3, PC4=PC4, clusters=clusters)
    datum <- na.omit(datum)
    model <- geeglm(cts ~ SEX + Genotype + PC1 + PC2 + PC3 + PC4, family="binomial",
                    data=datum, corstr="independence", id=clusters)
    coef.table <- coef(summary(model))
    wald <- as.numeric(coef.table$Wald[which(rownames(coef.table) == "Genotype")])
    return(c(SNP=as.character(snp.name),
             POS=as.character(snp.pos),
             beta=as.numeric(coef.table$Estimate[which(rownames(coef.table) == "Genotype")]),
             SE=as.numeric(coef.table$Std.err[which(rownames(coef.table) == "Genotype")]),
             OR=exp(as.numeric(coef.table$Estimate[which(rownames(coef.table) == "Genotype")])),
             Wald=wald,
             P.value=pchisq(wald,df=1,lower.tail=F)))  
  }
  
  CTS <- as.numeric(fam2[,6]) - 1
  FIDs <- as.numeric(fam2[,1])
    
  system.time(
  coefs <- foreach(i=1:num.snps, .combine=rbind) %dopar% {
      assoc.wrapper(Genotype=as.numeric(newped2[,i]), 
                  snp.name=as.character(map$V2[i]),
                  snp.pos=as.character(map$V4[i]),
                  pheno=CTS,
                  covariates=covar,
                  clusters=FIDs)
  }
  )
  print("Writing results to file")
  write.table(coefs,outname, quote=F, row.names=F, col.names=T)


