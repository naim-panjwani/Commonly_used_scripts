# Takes ped file with specified father and mother ID's and applies them to the fam file lacking these specifications

#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/CORRECTED_ANALYSIS")

ped <- read.table("17_REped_IIDs_updated.txt", header=T)
fam <- read.table("22_het_hap_snps_rm.fam",header=F)

assignParents <- function(iid) {
  # Given an IID, returns a vector of corresponding paternal and maternal IIDs
  pids <- as.character(ped[,3])
  mids <- as.character(ped[,4])
  index <- which(as.character(ped[,2]) == as.character(iid))
  return(c(pids[index],mids[index]))
}

temp <- sapply(as.character(fam[,2]), assignParents)
fam[,c(3,4)] <- cbind(sapply(temp, "[",1),
                      sapply(temp, "[",2))
fam[,6] <- 2

write.table(fam, "22_het_hap_snps_rm.fam", quote=F,row.names=F,col.names=F,sep=" ")