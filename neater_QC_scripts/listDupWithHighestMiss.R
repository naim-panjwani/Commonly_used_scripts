#!/usr/bin/R

library(doMC)
library(foreach)
registerDoMC(30)

args=(commandArgs(TRUE))

lmissfilename <- as.character(args[1])
outfilename <- as.character(args[2])

lmiss <- read.table(lmissfilename, header=T, stringsAsFactors=F)
dupnames <- gsub(".dup","", lmiss[grep("dup", lmiss$SNP),2])
dupnames2 <- c(dupnames, paste0(dupnames,".dup"))
lmiss2 <- subset(lmiss, lmiss$SNP %in% dupnames2)


dups_to_rm <- NULL
dups_to_rm <- foreach(i= 1:length(dupnames), .combine=c) %dopar% {
  indexes <- which(lmiss2[,2] %in% c(dupnames[i], paste0(dupnames[i],".dup")))
  tempmiss <- lmiss2[indexes,]
  return(tempmiss[which(tempmiss$N_MISS %in% max(tempmiss$N_MISS))[1],2])
}

write.table(dups_to_rm, outfilename, quote=F, row.names=F, col.names=F)

