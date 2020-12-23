# SLC26A9 plotting
# 3' to 5' (reverse gene)
# chr1:205,882,177-205,912,588

setwd("/Users/naim/Documents/Strug/R_code/")
genes <- read.table("slc26a9_refseq_genes_coordinates.txt",header=T, sep='\t')
dnase <- read.table("/Users/naim/Documents/Strug/SLC/dnase_26a9_extended.txt", header=T, sep='\t')
dnase <- dnase[,2:8]
chip <- read.table("/Users/naim/Documents/Strug/SLC/chip_v2_26a9.txt",header=T,sep='\t')

#============================================================================================
isInRange <- function(x, start, end) {
  # Given start and end numbers, returns TRUE if x is within start and end and FALSE otherwise
  result <- FALSE
  if(start > end) {
    temp <- start
    start <- end
    end <- temp
  }
  if(x >= start & x <= end) {
    result <- TRUE
  }
  return(result)
}
#-----------------------------------------------------------------
whichExons <- function(bp, exonStarts, exonEnds, reverse=FALSE) {
  # PRE: takes a range of bp positions (or just one position), exonStarts and exonEnds vectors of a gene
  # POST: returns a vector of the exons the range of bp positions fall into
  
  ifelse(length(bp)==1,bp <- c(bp,bp),
         ifelse(length(bp)!=2, stop("Not a valid base-pair range"), bp <- bp))
  if(length(exonStarts)!=length(exonEnds)) stop("Non-matching exon start and end coordinates")
  
  exons <- NULL
  for(i in 1:length(exonStarts)) {
    if(isInRange(bp[1],exonStarts[i],exonEnds[i]) | isInRange(bp[2],exonStarts[i],exonEnds[i]) |
         isInRange(exonStarts[i], bp[1],bp[2]) | isInRange(exonEnds[i],bp[1],bp[2])) {
      ifelse(!reverse, exons <- i, exons <- length(exonStarts)-i+1)
    }
  }
  return(exons)
}
#-----------------------------------------------------------------
whichIntrons <- function(bp, exonStarts, exonEnds, reverse=FALSE) {
  # PRE: takes a range of bp positions (or just one position), exonStarts and exonEnds vectors of a gene
  # POST: returns a vector of the introns the range of bp positions fall into
  ifelse(length(bp)==1,bp=c(bp,bp),
         ifelse(length(bp)!=2), stop("Not a valid base-pair range"), bp=bp)
  
}
#-----------------------------------------------------------------
colorExons <- function(bpStart,bpEnd, exonStarts, exonEnds, colour, reverse=FALSE) {
  # PRE: takes a range of bp positions and exon positions and desired colour 
  #      assumes a plot with exons is already drawn
  # POST: paints the exons in the bp range with color=colour
  
}
#============================================================================================
isoform1 <- genes[1,]
isoform2 <- genes[2,]
coord <- isoform1

upstreambp <- 5*1000 # 5kb
downstreambp <- 5*1000

genestart <- coord$txEnd + upstreambp
geneend <- coord$txStart - downstreambp
exonheight <- 0.2
exoncolor <- "blue"
UTRcolor <- "yellow"

exonStarts <- as.integer(unlist(strsplit(as.character(coord$exonStarts),",")))
exonEnds <- as.integer(unlist(strsplit(as.character(coord$exonEnds),",")))

png("/Users/naim/Desktop/26A9_plot.png", width=1440*2, height=1440, pointsize=24)
plot(x=c(geneend,genestart), y=c(0,0), xlab="chr1:205,882,177-205,912,588", type="l",
     main=paste("SLC26A9 Isoform 1\nRefSeq",as.character(coord$name)), ylab="")

# draw the exons
rect(exonStarts, -exonheight, exonEnds, exonheight, col=exoncolor)
UTR5start <- exonEnds[length(exonEnds)]
UTR5end <- coord$cdsEnd
UTR3start <- exonStarts[1]
UTR3end <- coord$cdsStart

# colorExons(UTR5end, UTR5start, UTRcolor)
# rect()


#rs3811424
points(205896948, 0.1, col="red", lwd=3)
#rs12133152
points(205895130, 0.1, col="darkred", lwd=3)
points(205894654, 0.1, col="pink", lwd=3)
points(205891942, 0.1, col="pink", lwd=3)
points(205891636, 0.1, col="pink", lwd=3)

#rs7512462
points(205899595, 0.1, col="darkgreen", lwd=3)
#rs1342063
points(205912859, 0.1, col="green", lwd=3)
points(205907872, 0.1, col="green", lwd=3)
points(205910080, 0.1, col="green", lwd=3)
points(205913073, 0.1, col="green", lwd=3)
points(205913848, 0.1, col="green", lwd=3)
points(205914757, 0.1, col="green", lwd=3)
points(205914885, 0.1, col="green", lwd=3)


# Draw the DNAseI peaks
dnase_maxheight <- 1
for(i in 1:nrow(dnase)) {
  lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=(c( ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)),
                                                 ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)) )),
        type='s',lwd=4)
#     lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=c(0.15,0.15),
#           type='s',lwd=4) 
}

legend("bottomright",legend = c("rs3811424", "rs7512462","Master DNAseI"), col=c("red","darkgreen","black"), pch=c(1,1,NA),lty=c(NA,NA,1),lwd=c(3,3,4))

dev.off()





# Zoom into rs12133152
png("/Users/naim/Desktop/26A9_plot_rs12133152.png", width=1440*2, height=1440, pointsize=24)
plot(x=c(205895130 - 5000, 205895130 + 5000), y=c(0,0), xlab="chr1:205895130 +/- 5kb", type="l",
     main=paste("SLC26A9 Isoform 1\nRefSeq",as.character(coord$name)), ylab="")
# draw the exons
rect(exonStarts, -exonheight, exonEnds, exonheight, col=exoncolor)
UTR5start <- exonEnds[length(exonEnds)]
UTR5end <- coord$cdsEnd
UTR3start <- exonStarts[1]
UTR3end <- coord$cdsStart

#rs3811424
points(205896948, 0.1, col="red", lwd=3)
#rs12133152
points(205895130, 0.1, col="darkred", lwd=3)
points(205894654, 0.1, col="pink", lwd=3)
points(205891942, 0.1, col="pink", lwd=3)
points(205891636, 0.1, col="pink", lwd=3)

#rs7512462
points(205899595, 0.1, col="darkgreen", lwd=3)
#rs1342063
points(205912859, 0.1, col="green", lwd=3)
points(205907872, 0.1, col="green", lwd=3)
points(205910080, 0.1, col="green", lwd=3)
points(205913073, 0.1, col="green", lwd=3)
points(205913848, 0.1, col="green", lwd=3)
points(205914757, 0.1, col="green", lwd=3)
points(205914885, 0.1, col="green", lwd=3)


# Draw the DNAseI peaks
dnase_maxheight <- 0.5
for(i in 1:nrow(dnase)) {
  lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=(c( ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)),
                                                          ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)) )),
        type='s',lwd=4)
  #     lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=c(0.15,0.15),
  #           type='s',lwd=4) 
}

legend("bottomright",legend = c("rs3811424", "rs12133152","rs7512462","Master DNAseI"), 
       col=c("red","darkred","darkgreen","black"), pch=c(1,1,1,NA),lty=c(NA,NA,NA,1),lwd=c(3,3,3,4))

dev.off()







# Zoom into rs7512462
png("/Users/naim/Desktop/26A9_plot_rs7512462.png", width=1440*2, height=1440, pointsize=32)
plot(x=c(205899595 - 1000, 205899595 + 20000), y=c(0,0), xlab="chr1:205899595 -1kb + 20kb", type="l",
     main=paste("SLC26A9 Isoform 1\nRefSeq",as.character(coord$name)), ylab="")
# draw the exons
rect(exonStarts, -exonheight, exonEnds, exonheight, col=exoncolor)
UTR5start <- exonEnds[length(exonEnds)]
UTR5end <- coord$cdsEnd
UTR3start <- exonStarts[1]
UTR3end <- coord$cdsStart

#rs7512462
points(205899595, 0.1, col="darkgreen", lwd=6)

points(205907872, 0.1, col="green", lwd=6)
points(205910080, 0.1, col="green", lwd=6)
points(205912859, 0.1, col="green", lwd=6)
points(205913073, 0.1, col="green", lwd=6)
points(205913848, 0.1, col="green", lwd=6)
points(205914757, 0.1, col="green", lwd=6)
points(205914885, 0.1, col="green", lwd=6)


# Draw the DNAseI peaks
dnase_maxheight <- 1.0
for(i in 1:nrow(dnase)) {
  lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=(c( ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)),
                                                          ((dnase$score[i] - min(dnase$score)) * dnase_maxheight) / (max(dnase$score) - min(dnase$score)) )),
        type='s',lwd=6)
  #     lines(x=c(dnase$chromStart[i],dnase$chromEnd[i]), y=c(0.15,0.15),
  #           type='s',lwd=4) 
}

legend("bottomright",legend = c("rs7512462","Master DNAseI"), 
       col=c("darkgreen","black"), pch=c(1,NA),lty=c(NA,1),lwd=c(6,6))

dev.off()

