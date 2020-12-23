# R script
args=(commandArgs(TRUE))

# Given an OLD bim file and text file with the updated SNP names, determines which ones
# differ and require changing. Use generate list_of_updated_snp_names",outname,".txt" file
# in conjuction with PLINK's --update-map <filename> --update-name

library(iterators)
library(parallel)
library(doParallel)
library(doMC)
library(foreach)
registerDoMC(30)
library(bigmemory)
library(biganalytics)
library(bigtabulate)

bimfilename <- as.character(args[1])
updated_snpnames_filename <- as.character(args[2])
outname <- as.character(args[3])

bim <- read.table(bimfilename,header=F)
bim <- cbind(bim, chr=bim[,1])
bim[bim[,7]==23,7]<-"X"
bim <- cbind(bim, pos=paste0("chr",as.character(bim[,7]),":",as.character(bim[,4])))

update <- read.table(updated_snpnames_filename,header=F)
update <- cbind(update,pos=paste0(as.character(update[,1]),":",as.character(update[,3])))

bim2 <- cbind(as.character(bim[,8]), as.character(bim[,2]),as.character(update[,4]))
colnames(bim2) <- c("POS","OLD","NEW")
bim_rows_unknown <- which(is.na(bim2[,3]))
bim2[bim_rows_unknown,3] <- bim2[bim_rows_unknown,2]
print("List of SNPs to be changed:")
bim_rows_to_change <- NULL
bim_rows_to_change <- foreach(i= 1:nrow(bim2), .combine=c) %dopar% {
  if(as.character(bim2[,2][i]) != as.character(bim2[,3][i])) {
    i
  }
}
(snps_to_change <- bim2[bim_rows_to_change,])
newbim <- bim
newbim[,2] <- bim2[,3]

write.table(snps_to_change, paste0("list_of_updated_snp_names",outname,".txt"), quote=F, row.names=F, col.names=T)
write.table(newbim[,1:6], outname, quote=F, row.names=F, col.names=F, sep="\t")
