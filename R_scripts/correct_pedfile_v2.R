# Takes ped file with specified father and mother ID's and applies them to the fam file lacking these specifications
# This version is immune to missing individuals within the corrected ped file (important in case of merged files)

#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/CORRECTED_ANALYSIS")

ped <- read.table("17_REped_IIDs_updated.txt", header=T)
fam <- read.table("27_merged_hapmap.fam",header=F)

assignParents <- function(iid) {
  # Given an IID, returns a vector of corresponding paternal and maternal IIDs
  result <- c("0","0")
  pids <- as.character(ped[,3])
  mids <- as.character(ped[,4])
  index <- which(as.character(ped[,2]) == as.character(iid))
  if(length(index)>0) result <- c(pids[index],mids[index])
  return(result)
}

temp <- sapply(as.character(fam[,2]), assignParents)
fam[,c(3,4)] <- cbind(temp[1,], temp[2,])
#fam[,6] <- 2

write.table(fam, "27_merged_hapmap.fam", quote=F,row.names=F,col.names=F,sep=" ")
