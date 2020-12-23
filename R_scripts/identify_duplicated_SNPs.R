#!/usr/bin/R
## Program will take duplicate positioned SNPs (not duplicate SNP names) and select the SNP with the lowest missing rate to keep ##
## Yes, this program will get rid of polymorphic SNPs as they may be in the same position
## PRE: .bim file from PLINK and .lmiss file for the exact same set of markers in .bim file
##      Enter the file names with the extensions please
## POST: returns the set of duplicated SNPs with the higher missing rate. If missing rates are the same, it tends to return the 2nd duplicated SNP in the list

args=(commandArgs(TRUE))
bimfile <- as.character(args[1])
lmissfile <- as.character(args[2])
out_prefix <- as.character(args[3]) # optional

if(is.na(out_prefix) | is.null(out_prefix)) out_prefix <- ""

library(data.table)

bim <- fread(bimfile, header=F, stringsAsFactors=F)
lmiss <- fread(lmissfile, header=T, stringsAsFactors=F)

MarkerName <- ifelse(bim$V1==0, bim$V2, paste0(bim$V1,":",bim$V4))
dupindexes <- c(which(duplicated(MarkerName)), which(duplicated(MarkerName, fromLast=T)))
dupindexes <- dupindexes[order(dupindexes)]

clusters <- data.frame(MarkerName=MarkerName[dupindexes], Index=dupindexes)
cluster <- 1
clusternum <- 1
marker <- clusters$MarkerName[1]
for(i in 2:nrow(clusters)) {
  currMarker <- clusters$MarkerName[i]
  if(currMarker == marker) {
    cluster <- c(cluster, clusternum)
  } else {
    marker <- clusters$MarkerName[i]
    clusternum <- clusternum + 1
    cluster <- c(cluster, clusternum)
  }
}

clusters <- cbind(clusters, cluster)


indexes_to_keep <- NULL
for(i in 1:max(clusters$cluster)) {
  clusti <- subset(clusters, clusters$cluster %in% i)
  currSurvivor <- clusti$Index[1]
  for(j in 2:nrow(clusti)) {
    if(lmiss$F_MISS[currSurvivor] > lmiss$F_MISS[clusti$Index[j]]) currSurvivor <- clusti$Index[j] 
  }
  indexes_to_keep <- c(indexes_to_keep, currSurvivor)
}

indexes_to_rm <- subset(dupindexes, !(dupindexes %in% indexes_to_keep))
snps_to_rm <- bim$V2[indexes_to_rm]

write.table(snps_to_rm, paste0(out_prefix, "-duplicated_SNPs_to_rm.txt"), quote=F, row.names=F, col.names=F)

