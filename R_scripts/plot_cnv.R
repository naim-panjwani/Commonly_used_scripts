#!/usr/bin/R
# Program to plot LRR and BAF similar to visualize_cnv.pl program
# Syntax: Rscript plot_cnv.R <filename>
# Example: Rscript plot_cnv.R sample.1054-306.chr13.44029637

# <filename> has the format with the following columns tab-delimited: Name, Chr, Position, sample.GType sample.B Allele Freq, sample.Log R Ratio
# The file should only have one region and be sorted in ascending order. Will plot everything into one plot


args=(commandArgs(TRUE))

filename <- as.character(args[1])

data <- read.table(filename, stringsAsFactors=F, header=T, sep="\t")

chrnum <- as.integer(data[1,2])
firstpos <- as.integer(data[1,3])
lastpos <- as.integer(data[nrow(data),3])
margin <- (lastpos - firstpos)*0.03
outname <- paste0(filename,".pdf")
pdf(outname)
par(mfrow=c(2,1))
plot(data[,3]/10^6, data[,6], main=filename, ylab="Log R Ratio", xlab=paste0("Chr", chrnum, " Position (Mb)"),
      xlim=c((firstpos - margin)/10^6, (lastpos + margin)/10^6), ylim=c(-1.0, 1.0),
      pch=20,cex=0.3)
abline(h=0,col="red",lty=2,lwd=0.3)
plot(data[,3]/10^6, data[,5], main=filename, ylab="B Allele Freq", xlab=paste0("Chr", chrnum, " Position (Mb)"),
      xlim=c((firstpos - margin)/10^6, (lastpos + margin)/10^6), ylim=c(0.0, 1.0),
      pch=20,cex=0.3)
abline(h=0.5,col="red",lty=2,lwd=0.3)
dev.off()

bitmap(file=paste0(filename,".jpg"), type="jpeg", width=3, height=3, res=300, pointsize=6)
par(mfrow=c(2,1))
plot(data[,3]/10^6, data[,6], main=filename, ylab="Log R Ratio", xlab=paste0("Chr", chrnum, " Position (Mb)"),
      xlim=c((firstpos - margin)/10^6, (lastpos + margin)/10^6), ylim=c(-1.0, 1.0),
      pch=20,cex=0.3)
abline(h=0,col="red",lty=2,lwd=0.3)
plot(data[,3]/10^6, data[,5], main=filename, ylab="B Allele Freq", xlab=paste0("Chr", chrnum, " Position (Mb)"),
      xlim=c((firstpos - margin)/10^6, (lastpos + margin)/10^6), ylim=c(0.0, 1.0),
      pch=20,cex=0.3)
abline(h=0.5,col="red",lty=2,lwd=0.3)
dev.off()


