#!/usr/bin/R
args=(commandArgs(TRUE))

# Input: chr, start, end, bimfilename
# Output: # of SNPs in the given range

chr=as.character(args[1])
if(chr == "X") {
  chr <- as.numeric(23)
} else {
  chr <- as.numeric(args[1])
}
start=as.numeric(args[2])
end=as.numeric(args[3])
bimfilename=as.character(args[4])

bimfile <- read.table(bimfilename, header=F)

num_snps <- 0
bim_subset <- subset(bimfile, bimfile[,1] == chr & bimfile[,4] >= start & bimfile[,4] <= end)
num_snps <- nrow(bim_subset)

write.table(num_snps,"tmp_num_snps.txt",quote=F,row.names=F,col.names=F)

