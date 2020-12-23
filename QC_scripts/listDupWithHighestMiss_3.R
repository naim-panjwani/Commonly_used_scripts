#!/usr/bin/R

# PRE: ensure all SNP names have unique names (use add-dup.R prior to running this)
# POST: identifies best SNP among duplicates and creates a list of worst duplicates to exclude according to call rate; 
#  use PLINK to exclude the set output

# Note: a duplicate SNP is one that has same chr:pos:ref:alt combination

args=(commandArgs(TRUE))

bimfilename <- as.character(args[1])
lmissfilename <- as.character(args[2])
listpolymorphic <- as.character(args[3]) # True of false?
outfilename <- as.character(args[4])


library(data.table)

bim <- fread(bimfilename, header=F, stringsAsFactors=F)
lmiss <- fread(lmissfilename, header=T, stringsAsFactors=F)

listpolymorphic <- ifelse(as.character(listpolymorphic) %in% c("FALSE", "False", "F", "0", "N", "No", "NO"), FALSE, TRUE) 

flipAllele <- function(a) {
  newa <- ifelse(a=="A", "T", ifelse(a=="T", "A", ifelse(a=="G","C", ifelse(a=="C", "G", a))))
  return(newa)
}

sameAlleles <- function(a1, a2) {
  if(length(unique(a1))==1 & length(unique(a2))==1) return(TRUE)
  firstpair <- c(a1[1], a2[1])
  for(i in 2:length(a1)) {
    pairi <- c(a1[i], a2[i])
    if(length(subset(pairi, pairi %in% firstpair)) != 2) {
      pairi <- c(flipAllele(pairi[1]), flipAllele(pairi[2]))
      if(length(subset(pairi, pairi %in% firstpair)) != 2) return(FALSE)
    }
  }
  return(TRUE)
}

getUniqueAllelePairs <- function(a1, a2) {
  firstpair <- c(a1[1], a2[1])
  uniquepairs <- data.table(a1=firstpair[1], a2=firstpair[2])
  for(i in 2:length(a1)) {
    pairi <- c(a1[i], a2[i])
    if(!(sameAlleles( a1=c(pairi[1], uniquepairs[,a1]), a2=c(pairi[2], uniquepairs[,a2]) ) )) uniquepairs <- rbind(uniquepairs, data.table(a1[i], a2[i]), use.names=F)
  }
  return(uniquepairs)
}


altsnpname <- sapply(1, function(i) paste0(bim[,V1], ":", bim[,V4], ":", bim[,V5], ":", bim[,V6]))[,1]
positions <- sapply(1, function(i) paste0(bim[,V1],":", bim[,V4]))[,1]
df <- cbind(bim[,.(V1,V2,V4,V5,V6)], lmiss[,.(N_MISS,N_GENO,F_MISS)], positions, altsnpname)

dups <- unique(altsnpname[duplicated(altsnpname)])
duppos <- unique(positions[duplicated(positions)])
excl <- NULL
counter <- 1
for(dupos in duppos) {
  df_subset <- subset(df, df[,positions] %in% dupos)
  if(sameAlleles(df_subset[,V5], df_subset[,V6])) {
    keeper_i <- which(df_subset$N_MISS %in% min(df_subset$N_MISS))[1]
    nonkeepers <- df_subset[-keeper_i, V2]
    excl <- c(excl, nonkeepers)
  }
  else if (listpolymorphic) {
    excl <- c(excl, df_subset[,V2])
  }
  else {
    print(paste("entered here for",dupos))
    uniqueAllelePairs <- getUniqueAllelePairs(df_subset[,V5], df_subset[,V6])
    for(i in 1:nrow(uniqueAllelePairs)) {
      currAllelePair <- c(uniqueAllelePairs[i,a1], uniqueAllelePairs[i,a2])
      df_subset2 <- NULL
      for(j in 1:nrow(df_subset)) {
        allelestmp <- data.table(a1=currAllelePair[1], a2=currAllelePair[2])
        allelestmp <- rbind(allelestmp, data.table(a1=df_subset[j, V5], a2=df_subset[j, V6]), use.names=F)
        if(sameAlleles(allelestmp[,a1], allelestmp[,a2])) df_subset2 <- rbind(df_subset2, df_subset[j,])
      }
      if(sameAlleles(df_subset2[,V5], df_subset2[,V6])) {
        keeper_i <- which(df_subset2$N_MISS %in% min(df_subset2$N_MISS))[1]
        nonkeepers <- df_subset2[-keeper_i, V2]
        if(length(nonkeepers)>0) excl <- c(excl, nonkeepers)
      }
    }
  }
  print(paste("finished",counter,"of",length(duppos)))
  counter <- counter+1
}

if(is.null(excl)) print("No duplicate SNPs were found")

write.table(excl, outfilename, quote=F, row.names=F, col.names=F)

