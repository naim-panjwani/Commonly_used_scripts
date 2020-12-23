# R script to find and list all duplicate SNPs
# You may want to run it twice or thrice after the 1st round of removal

args=(commandArgs(TRUE))
#setwd("/Users/naim/Documents/Strug/141202-Spit_for_Science_all_samples/QC/")

bimfilename <- as.character(args[1])
lmissfilename <- as.character(args[2])
outname <- as.character(args[3])

bim <- read.table(bimfilename,header=F)
lmiss <- read.table(lmissfilename, header=T)

dup.snps <- NULL #Will store all duplicated SNPs
worst.dup.snp <- NULL # Will store the list with lowest call rate for each duplicated SNP
for(chr in 1:23) {
  snps <- subset(bim, bim[,1] == chr)
  num.snps <- nrow(snps)
  indices <- which(duplicated(snps[,4], fromLast=TRUE))
  indices2 <- which(duplicated(snps[,4]))
  all_indices <- NULL
  if(length(indices)>0) {
    for(i in 1:length(indices)) all_indices <- c(all_indices, indices[i],indices2[i])
  }
  dup.snps <- rbind(dup.snps,snps[all_indices,])
  
  miss <- subset(lmiss, lmiss[,1] == chr)
  temp <- NULL
  if(length(indices)>0) {
    temp <- cbind(as.numeric(miss[indices,5]), as.numeric(miss[indices2,5]), indices, indices2)
    itemp <- 3+(temp[,1]<temp[,2])
  }
  worst.dup.snp.indices <- NULL
  if(!is.null(temp)) {
    for(i in 1:nrow(temp)) worst.dup.snp.indices <- c(worst.dup.snp.indices, as.numeric(temp[i,itemp[i]]))
    worst.dup.snp <- rbind(worst.dup.snp, miss[worst.dup.snp.indices,])
  }
}
write.table(dup.snps,paste0(outname,"_duplicated_snps_details.txt"), quote=F, row.names=F, col.names=F)
write.table(worst.dup.snp, paste0(outname,"_worst_duplicated_snps.txt"), quote=F, row.names=F, col.names=F)

