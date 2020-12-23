plot_locuszoom <- function(chr, start, end, ld_filename, lead_snpname, transparency, superimpose=TRUE, yaxis="left") {
  GWAS_pvalues <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% chr & 
                           All_GWAS_pvalues$pos >= start & All_GWAS_pvalues$pos <= end)
  ld <- read.table(ld_filename,
                   header=T, stringsAsFactors = F)
  names(ld)[which(names(ld) %in% "SNP_B")] <- "MarkerName"
  ld_mergedset <- merge(GWAS_pvalues, ld[,c(6,7)], by="MarkerName")
  
  lead_snp <- subset(ld_mergedset, ld_mergedset$MarkerName %in% lead_snpname) #rs7549173
  r2_80_group <- subset(ld_mergedset, ld_mergedset$R2>=0.8)
  r2_60_group <- subset(ld_mergedset, ld_mergedset$R2>=0.6 & ld_mergedset$R2<0.8)
  r2_40_group <- subset(ld_mergedset, ld_mergedset$R2>=0.4 & ld_mergedset$R2<0.6)
  r2_20_group <- subset(ld_mergedset, ld_mergedset$R2>=0.2 & ld_mergedset$R2<0.4)
  r2_lt_20_group <- subset(ld_mergedset, ld_mergedset$R2<0.2)
  no_ld_info_snps <- subset(GWAS_pvalues, !(GWAS_pvalues$MarkerName %in% ld_mergedset$MarkerName))
  
  if(superimpose==T) {
    par(new=T)
    plot(c(start/10^6, end/10^6), c(0,-log10(lead_snp$p)+1), axes=F, xlab=NA, ylab=NA, type="n")
    if(yaxis %in% "right") {
      axis(side=4)
      mtext("GWAS -log10(P-values)", side=4, line=2)
    }
    text(x=lead_snp$pos/10^6,y=-log10(lead_snp$p)+0.5, labels=lead_snp$snp)
    points(r2_80_group$pos/10^6, -log10(r2_80_group$p), pch=16, col=alpha("red",transparency),cex=1.5)
    points(lead_snp$pos/10^6, -log10(lead_snp$p), pch=23, col=alpha("purple",transparency), bg=alpha("purple",transparency),cex=1.5)
    points(r2_60_group$pos/10^6, -log10(r2_60_group$p), pch=16, col=alpha("orange",transparency),cex=1.5)
    points(r2_40_group$pos/10^6, -log10(r2_40_group$p), pch=16, col=alpha("green",transparency),cex=1.5)
    points(r2_20_group$pos/10^6, -log10(r2_20_group$p), pch=16, col=alpha("cyan",transparency),cex=1.5)
    points(r2_lt_20_group$pos/10^6, -log10(r2_lt_20_group$p), pch=16, col=alpha("blue",transparency),cex=1.5)
    points(no_ld_info_snps$pos/10^6, -log10(no_ld_info_snps$p), pch=16, col=alpha("black",transparency),cex=1.5)
  } else {
    plot(x=c(plot_region$start/10^6, plot_region$end/10^6), y=c(0,-log10(lead_snp$p)+1), 
         xlab=paste0("chr",chr," coordinates (Mb)"), ylab="GWAS -log10(P-values)",type="n")
    if(yaxis %in% "right") {
      axis(side=4)
      mtext("GWAS -log10(P-values)", side=4, line=2)
    }
    text(x=lead_snp$pos/10^6,y=-log10(lead_snp$p)+0.5, labels=lead_snp$snp)
    points(r2_80_group$pos/10^6, -log10(r2_80_group$p), pch=16, col=alpha("red",transparency),cex=1.5)
    points(lead_snp$pos/10^6, -log10(lead_snp$p), pch=23, col=alpha("purple",transparency), bg=alpha("purple",transparency),cex=1.5)
    points(r2_60_group$pos/10^6, -log10(r2_60_group$p), pch=16, col=alpha("orange",transparency),cex=1.5)
    points(r2_40_group$pos/10^6, -log10(r2_40_group$p), pch=16, col=alpha("green",transparency),cex=1.5)
    points(r2_20_group$pos/10^6, -log10(r2_20_group$p), pch=16, col=alpha("cyan",transparency),cex=1.5)
    points(r2_lt_20_group$pos/10^6, -log10(r2_lt_20_group$p), pch=16, col=alpha("blue",transparency),cex=1.5)
    points(no_ld_info_snps$pos/10^6, -log10(no_ld_info_snps$p), pch=16, col=alpha("black",transparency),cex=1.5)
  }
}



# Example:
  # ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7549173_slc26a9.ld"
  # lead_snpname <- "chr1:205906897"
  # transparency <- 0.2
  # plot_locuszoom("1",205860000,205940000, ld_filename, lead_snpname, transparency, superimpose=F, yaxis="left")


