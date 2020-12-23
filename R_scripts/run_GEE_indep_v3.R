#!/hpf/tools/centos6/R/3.1.1/bin/R
# GEE Independence
# CTS = Genotype + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7
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
registerDoMC(50)
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
  mpermute(ped, order=as.numeric(newfam[,1]))
  newped2 <- sub.big.matrix(ped, firstRow = 1, lastRow = nrow(ped), firstCol = 7, lastCol = ncol(ped))
   print("Data re-ordering done")
  
  
  # 3. Generate covariates variable
  print("Generating covariates")
  covar <- NULL
  covar <- cbind(fam2[,1:2]) # FID/IID
  sex <- fam2[,5]
  PCs <- newpca2[,7:13] # 7 PCs 
  covar <- cbind(covar, sex, PCs)
  names(covar) <- c("FID", "IID", "SEX", "PC1", "PC2","PC3", "PC4", "PC5", "PC6", "PC7")

  snpnames <- as.character(map$V2)
  chr <- map[,1]
  CTS <- as.numeric(fam2[,6]) - 1
  FIDs <- as.numeric(fam2[,1])

  # 5. Model the association
  
  print("Running GEE independence association model")
  assoc.wrapper <- function(Genotype, snp.chr, snp.pos, snp.name, pheno, covariates, clusters) {
    SEX <- as.numeric(covar$SEX) -1
    PC1 <- covariates$PC1
    PC2 <- covariates$PC2
    PC3 <- covariates$PC3
    PC4 <- covariates$PC4
    PC5 <- covariates$PC5
    PC6 <- covariates$PC6
    PC7 <- covariates$PC7
    datum <- data.frame(cts=pheno, SEX=SEX, Genotype=Genotype, PC1=PC1, PC2=PC2, PC3=PC3, PC4=PC4, PC5=PC5, PC6=PC6, PC7=PC7, clusters=clusters)
    datum <- na.omit(datum)
    tryCatch(
    {
      model <- geeglm(cts ~ SEX + Genotype + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7, family="binomial",
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
    }
    , error = function(e) {return(c(CHR=as.character(snp.chr), POS=as.character(snp.pos), SNP=as.character(snp.name), beta=NA, SE=NA, OR=NA, Wald=NA, P.value=NA))}
   ) 
  }
  
    
  system.time(
  coefs <- foreach(i=1:num.snps, .combine=rbind) %dopar% {
      assoc.wrapper(Genotype=as.numeric(newped2[,i]),
                    snp.chr=as.numeric(map[i,1]),
                    snp.pos=as.character(map$V4[i]),
                    snp.name=as.character(map$V2[i]),
                    pheno=CTS,
                    covariates=covar,
                    clusters=FIDs)
  }
  )
  print("Writing results to file")
  write.table(coefs,outname, quote=F, row.names=F, col.names=T)


