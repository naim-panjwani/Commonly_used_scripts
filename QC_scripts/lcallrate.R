#!/usr/bin/R
# Plots per-SNP call rate vs percentile for all SNPs and those above call rate > 90% (prefix_ + lcallrate.pdf)
# SNPs with low call rate (<=90%) are written to file (prefix_ + lowCallSNPs.txt)
# Syntax: Rscript lcallrate.R <filename.lmiss> <prefix (optional)>
# Outputs: (prefix_)lcallrate.png, (prefix_)lowCallSNPs.txt

args=(commandArgs(TRUE))

library(data.table)

lmissfilename <- as.character(args[1])
prefix <- as.character(args[2])
if(is.na(prefix)) {
  prefix <- ""
} else {
  prefix <- paste0(prefix,"_")
}
lmiss <- fread(lmissfilename, header=TRUE, stringsAsFactors=F)
lmiss <- na.omit(lmiss)

#pdf(paste0(prefix ,"lcallrate.pdf"))
png(paste0(prefix ,"lcallrate.png"))
par(mfrow=c(1,2))
callrate <- (1-lmiss$F_MISS)*100
quants <- as.numeric(gsub("%","",names(quantile(callrate, probs=seq(0,1,1/(length(callrate)-1))))))
callrateSorted <- quantile(callrate, probs=seq(0,1,1/(length(callrate)-1)))
plot(quants, callrateSorted, xlab="Percentile (%)", ylab="SNP Call Rate (%)", main="All SNPs", type="n")
points(subset(quants, callrateSorted<=90), subset(callrateSorted, callrateSorted<=90), col="red")
points(subset(quants, callrateSorted>90), subset(callrateSorted, callrateSorted>90), col="black")

callrate <- (1-subset(lmiss$F_MISS, lmiss$F_MISS<=0.1))*100
quants <- as.numeric(gsub("%","",names(quantile(callrate, probs=seq(0,1,1/(length(callrate)-1))))))
callrateSorted <- quantile(callrate, probs=seq(0,1,1/(length(callrate)-1)))
plot(quants, callrateSorted, xlab="Percentile (%)", ylab="SNP Call Rate (%)", main="SNPs with call rate > 90%", type="n")
points(subset(quants, callrateSorted<=90), subset(callrateSorted, callrateSorted<=90), col="red")
points(subset(quants, callrateSorted>90), subset(callrateSorted, callrateSorted>90), col="black")
dev.off()

highmissSNPs <- subset(lmiss, lmiss$F_MISS>=0.1)
write.table(highmissSNPs, paste0(prefix, "lowCallSNPs.txt"), quote=F, row.names=F, col.names=T)


