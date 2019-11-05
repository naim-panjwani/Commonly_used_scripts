#!/usr/bin/R
# Plots per-sample call rate vs percentile for all individuals and those above call rate > 90% (prefix_ + icallrate.pdf)
# Individuals with low call rate (<=90%) are written to file (prefix_ + lowCallInd.txt)
# Syntax: Rscript icallrate.R <filename.imiss> <prefix (optional)>
# Outputs: (prefix_)icallrate.pdf, (prefix_)lowCallInd.txt

args=(commandArgs(TRUE))

imissfilename <- as.character(args[1])
prefix <- as.character(args[2])
if(is.na(prefix)) {
  prefix <- ""
} else {
  prefix <- paste0(prefix,"_")
}
imiss <- read.table(imissfilename, header=TRUE, stringsAsFactors=F)

pdf(paste0(prefix ,"icallrate.pdf"))
par(mfrow=c(1,2))
callrate <- (1-imiss$F_MISS)*100
quants <- as.numeric(gsub("%","",names(quantile(callrate, probs=seq(0,1,1/(length(callrate)-1))))))
callrateSorted <- quantile(callrate, probs=seq(0,1,1/(length(callrate)-1)))
plot(quants, callrateSorted, xlab="Percentile (%)", ylab="Individual Call Rate (%)", main="All individuals", type="n")
points(subset(quants, callrateSorted<=90), subset(callrateSorted, callrateSorted<=90), col="red")
points(subset(quants, callrateSorted>90), subset(callrateSorted, callrateSorted>90), col="black")

callrate <- (1-subset(imiss$F_MISS, imiss$F_MISS<=0.1))*100
quants <- as.numeric(gsub("%","",names(quantile(callrate, probs=seq(0,1,1/(length(callrate)-1))))))
callrateSorted <- quantile(callrate, probs=seq(0,1,1/(length(callrate)-1)))
plot(quants, callrateSorted, xlab="Percentile (%)", ylab="Individual Call Rate (%)", main="Individuals with call rate > 90%", type="n")
points(subset(quants, callrateSorted<=90), subset(callrateSorted, callrateSorted<=90), col="red")
points(subset(quants, callrateSorted>90), subset(callrateSorted, callrateSorted>90), col="black")
dev.off()

highmissInd <- subset(imiss, imiss$F_MISS>=0.1)
flagmissInd <- subset(imiss, imiss$F_MISS>=0.01 & imiss$F_MISS<0.1)
write.table(highmissInd, paste0(prefix, "lowCallInd.txt"), quote=F, row.names=F, col.names=T)
write.table(flagmissInd, paste0(prefix, "flagged_90-99callRateInds.txt"), quote=F, row.names=F, col.names=T)

