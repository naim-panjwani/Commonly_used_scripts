cts<-read.table("CTS_cases.bim",header=F)
ctrls<-read.table("chr11.block4.11p13.dosage",header=F)[,1:10]
common_snps<-NULL
for(i in 1:dim(cts)[1]) { if(as.character(cts$V2[i]) %in% as.character(ctrls$V2)) common_snps<-c(common_snps,as.character(cts$V2)[i]) }
write.table(common_snps,"common_snps.txt",quote=F,row.names=F,col.names=F)
