#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/141020-GWAS_GEE_TOR_CF_Ctrls_CTS_exome")

columns_to_extract <- read.table("TOR_CF_600ind_details.txt") #V4 is column of interest
colnames(columns_to_extract) <- c("FID","IID","Sex","dosage_columns")

# Write down how many blocks per chromosome folder
# eg. chr1 has 25 blocks
block_lengths <- c(25, 25, 20, 20, 19, 
                   18, 16, 15, 15, 14,
                   14, 14, 12, 11, 11,
                   10, 9, 8, 6, 7,
                   5, 6, 16)
block_starts <- c(1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,2,3,3,
                  1,1,1,1,1,
                  2,2,1)

dosageTOgenotype<- function(minorAlleleCount, minorAllele, majorAllele) {
  # PRE: minorAlleleCount is 0,1,2
  #      minorAllele is the minorAllele nucleotide: A,C,G,T
  #      majorAllele is the majorAllele nucleotide: A,C,G,T
  # POST: gives back corresponding genotype
  # Examples: dosageTOgenotype(0,"A","T") returns TT
  
  geno <- NULL
  geno<-ifelse(minorAlleleCount==0, paste(majorAllele,majorAllele,sep=" "),
               ifelse(minorAlleleCount==1, paste(minorAllele,majorAllele,sep=" "),
                      ifelse(minorAlleleCount==2,paste(minorAllele,minorAllele,sep=" "),NA)))
  return(geno)
}


i = 1
#for (i in 1:23) { # for each chromosome
  
  print(paste("chr",i))
  for (j in block_starts[i]:block_lengths[i]) { # for each block
    setwd(paste("./GWAS1_imputation/chr",i,sep=""))
    print(paste("Reading chr",i,"block",j))
    dosage.ctrls <- read.table(paste("chr",i,".block",j,".dosage.gz",sep=""))[,c(1:5,columns_to_extract$dosage_columns)]
    colnames(dosage.ctrls) <- c("CHR","SNP","BP","A1","A2", as.character(columns_to_extract$IID))
    dosage_geno_data <- round(dosage.ctrls[,6:dim(dosage.ctrls)[2]])
    
    print(paste("Converting genotype data chr",i,"block",j))
    new_dosage<-matrix(nrow=nrow(dosage_geno_data),ncol=ncol(dosage_geno_data))
    for(k in 1:nrow(dosage_geno_data)) {
      minorAllele<-as.character(dosage.ctrls$A1[k])
      majorAllele<-as.character(dosage.ctrls$A2[k])
      new_dosage[k,]<-sapply(dosage_geno_data[k,],dosageTOgenotype,minorAllele,majorAllele)
    }
    colnames(new_dosage)<-as.character(columns_to_extract$IID)
    
    new_dosage.ctrls.tped<-cbind(dosage.ctrls[,1:2],cM=rep(0,nrow(dosage.ctrls)),dosage.ctrls[,3],new_dosage)
    tfam<-cbind(FID=as.character(columns_to_extract$FID),IID=as.character(columns_to_extract$IID),mid=0,pid=0,sex=columns_to_extract$SEX,phenotype=1)
    
    setwd("/home/naim/UKRE/UKRE_assoc/141020-GWAS_GEE_TOR_CF_Ctrls_CTS_exome")
    filename <- paste("chr",i,"block",j,sep="")
    print(paste("Saving chr",i,"block",j))
    write.table(new_dosage.ctrls.tped, file = paste("02_",filename,".tped",sep=""),quote=F,row.names=F,col.names=F)
    write.table(tfam, file = paste("02_",filename,".tfam",sep=""),quote=F,row.names=F,col.names=F)    
  }
#}



