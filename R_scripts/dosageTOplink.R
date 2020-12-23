setwd("/Users/naim panjwani/Documents/Strug/UKRE/UKRE_assoc/TOR_CF_Ctrls")

pedfile<-read.table("CTS_cases_common_snps.ped",header=F)
mapfile<-read.table("CTS_cases_common_snps.map",header=F)
dosage.ctrls<-read.table("chr11.block4.11p13.dosage",header=F)
pheno<-read.csv("phenotype.gwas1.csv",header=T)
dosage.ctrls<-subset(dosage.ctrls, as.character(dosage.ctrls$V2) %in% as.character(mapfile$V2))
num.snps<-dim(mapfile)[1]
num.individuals<-dim(pedfile)[1]
snp.names <- mapfile$V2

dosage_data<-round(dosage.ctrls[,6:dim(dosage.ctrls)[2]])
colnames(dosage_data)<-as.character(pheno$IID)
colnames(dosage.ctrls)<-c("CHR","SNP","BP","A1","A2",as.character(pheno$IID))

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

new_dosage<-matrix(nrow=nrow(dosage_data),ncol=ncol(dosage_data))
for(i in 1:nrow(dosage_data)) {
  minorAllele<-as.character(dosage.ctrls$A1[i])
  majorAllele<-as.character(dosage.ctrls$A2[i])
  new_dosage[i,]<-sapply(dosage_data[i,],dosageTOgenotype,minorAllele,majorAllele)
}
colnames(new_dosage)<-as.character(pheno$IID)

new_dosage.ctrls.tped<-cbind(dosage.ctrls[,1:2],cM=rep(0,nrow(dosage.ctrls)),dosage.ctrls[,3],new_dosage)
tfam<-cbind(FID=as.character(pheno$FID),IID=as.character(pheno$IID),mid=0,pid=0,sex=pheno$SEX,phenotype=1)

  write.table(new_dosage.ctrls.tped,"01_TOR_CF_Ctrls.tped",quote=F,row.names=F,col.names=F)
  write.table(tfam, "01_TOR_CF_Ctrls.tfam",quote=F,row.names=F,col.names=F)
  
  tor_ind<-subset(tfam,as.character(pheno$SITE) %in% "TOR")[,1:2]
  write.table(tor_ind,"TOR_individuals.txt",quote=F,row.names=F,col.names=F)
  
  tor_bc_ind <- subset(tfam, as.character(pheno$COHORT) %in% "TOR-BC")
  write.table(tor_bc_ind,"TOR_BC_individuals.txt", quote=F,row.names=F,col.names=F)
  
