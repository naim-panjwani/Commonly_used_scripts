setwd("/Users/naim/Documents/Strug/CF/MI_GWAS_tissue_enrichment/")
#gtex_dir <- "/Users/naim/Documents/Strug/GTEx_data/original/subsets/"
gtex_dir <- "/Users/naim/Documents/Strug/GTEx_data/eQTL_plots_v7_update/"

#========================================= Load Data ======================================================
library(scales)
par(oma=c(0,0,0,4))
# par(mar=c(5,4,4,4)+0.1)


###### CHR7 is outdated for MI_FDR_10-5_SNPs.txt FILE #####
# MI_suggestive_pvalues <- read.table("MI_FDR_10-5_SNPs.txt", header=F, stringsAsFactors = F)
# MI_suggestive_pvalues$V4 <- ifelse(MI_suggestive_pvalues$V4==0,1e-16,MI_suggestive_pvalues$V4)
# MI_suggestive_pvalues <- cbind(MarkerName=paste0(ifelse(MI_suggestive_pvalues$V1==23,"X",MI_suggestive_pvalues$V1),
#                                                  ":",MI_suggestive_pvalues$V3), MI_suggestive_pvalues)

# MI overall
# All_GWAS_pvalues <- read.table("FDR_noNA_G12.txt", header=T, stringsAsFactors = F)  ### This is the older data imputed with 1KG only
#All_GWAS_pvalues <- read.table("FDR_out_G12_hybrid_ref_imputed_assoc.txt", header=T, stringsAsFactors = F) ### This is the association redone with the 1KG+101CF hybrid imputation
All_GWAS_pvalues <- read.table("hybrid_assoc_regions_of_interest.txt", header=T, stringsAsFactors = F) ### This is the association redone with the 1KG+101CF hybrid imputation (corrected round)
All_GWAS_pvalues <- cbind.data.frame(All_GWAS_pvalues, MarkerName=paste0("chr",All_GWAS_pvalues$chr,":",All_GWAS_pvalues$pos), stringsAsFactors=F)
if(length(which(All_GWAS_pvalues$p==0))!=0) All_GWAS_pvalues[which(All_GWAS_pvalues$p==0),]$p <- 2.2e-16 
All_GWAS_pvalues <- na.omit(All_GWAS_pvalues)

#### No need to update chr7 anymore as the whole association file is updated above
# update chr7:
# chr7 <- read.table("newChr7MI.txt", header=T, stringsAsFactors = F)
# colnames(chr7) <- colnames(All_GWAS_pvalues)
# chr7 <- na.omit(chr7)
# All_GWAS_pvalues <- rbind(subset(All_GWAS_pvalues, !(All_GWAS_pvalues$chr %in% 7)), chr7)

# Males lung (SLC6A14)
lung_assoc_pvalues <- read.table("chrX_lung_association_partad_femaleADD.txt", header=T, stringsAsFactors = F)
lung_assoc_pvalues <- lung_assoc_pvalues[,c('SNP','CHR','BP','P_value')] # Male 0/2 coding association additive model
names(lung_assoc_pvalues) <- c("snp","chr","pos","p")
lung_assoc_pvalues <- cbind(lung_assoc_pvalues, MarkerName=paste0("chr",lung_assoc_pvalues$chr,":",lung_assoc_pvalues$pos))
lung_assoc_pvalues <- na.omit(lung_assoc_pvalues)

# Females lung (SLC6A14)
lung_assoc_pvalues_females <- read.table("chrX_lung_association_partad_femaleADD.txt", header=T, stringsAsFactors = F)
lung_assoc_pvalues_females <- lung_assoc_pvalues_females[,c('SNP','CHR','BP','P_value2')]
names(lung_assoc_pvalues_females) <- c("snp","chr","pos","p")
lung_assoc_pvalues_females <- cbind(lung_assoc_pvalues_females, MarkerName=paste0("chr", lung_assoc_pvalues_females$chr,":",lung_assoc_pvalues_females$pos))
lung_assoc_pvalues_females <- na.omit(lung_assoc_pvalues_females)

# Males+females lung (SLC6A14) -- these results are from our data -- guessing this is done by Bowei (additive model, mega-GWAS)
# lung_assoc_overall <- read.table("lungAssociation_6A14region_bothSex.dsg", header=T, stringsAsFactors = F)
# names(lung_assoc_overall) <- c("chr","pos","snp","Beta","SE","p","AF","n")
# lung_assoc_overall <- lung_assoc_overall[,c('snp','chr','pos','p')]
# lung_assoc_overall <- cbind(lung_assoc_overall, MarkerName=paste0("chr",lung_assoc_overall$chr,":",lung_assoc_overall$pos))
# lung_assoc_overall <- na.omit(lung_assoc_overall)


# Using the published GWAS instead (by UNC):
lung_assoc_overall_UNC <- read.table("gwas.public_lung_assoc_6a14.txt", header=T, stringsAsFactors = F)
lung_assoc_overall_UNC <- lung_assoc_overall_UNC[,c(3,1,2,7)]
names(lung_assoc_overall_UNC) <- c("snp","chr","pos","p")
lung_assoc_overall_UNC <- cbind(lung_assoc_overall_UNC, MarkerName=paste0("chr",lung_assoc_overall_UNC$chr,":",lung_assoc_overall_UNC$pos))

#lung_assoc_overall <- lung_assoc_overall_UNC


# MI Males (SLC6A14)
#MI_assoc_chrX <- read.table("MI_assoc_chrX_by_sex_Additive.txt", header=T, stringsAsFactors = F) # needs updating from correct hybrid imputation
MI_assoc_males <- MI_assoc_chrX[,c('SNP','Chr','BP','P_Male')]
names(MI_assoc_males) <- c("snp","chr","pos","p")
MI_assoc_males <- cbind(MI_assoc_males, MarkerName=paste0("chr",MI_assoc_males$chr,":",MI_assoc_males$pos))
MI_assoc_males <- na.omit(MI_assoc_males)

# MI Females (SLC6A14)
MI_assoc_females <- MI_assoc_chrX[,c('SNP','Chr','BP','P_Female')]
names(MI_assoc_females) <- c("snp","chr","pos","p")
MI_assoc_females <- cbind(MI_assoc_females, MarkerName=paste0("chr",MI_assoc_females$chr,":",MI_assoc_females$pos))
MI_assoc_females <- na.omit(MI_assoc_females)

#### Enhancer and promoter states from Broad Institute: -log10P DNAse > 2 and using the 15-state chromHMM ####
#relevant_tissues <- read.table("per_locus_contrast_tests/tissue_contrasts/relevant_tissues_list.txt", sep="\t", stringsAsFactors = F) ## too many; using 2nd list
relevant_tissues <- read.table("per_locus_contrast_tests/tissue_contrasts/relevant_tissues_list3.txt", sep="\t", stringsAsFactors = F) ## less tracks than v1

load("/Users/naim/Documents/chromHMM_state_calls_DNaselog10_gt_2_enhancers_15state.rdata")
enhancers <- max_states
load("/Users/naim/Documents/chromHMM_state_calls_DNaselog10_gt_2_promoters_15state.RData")
promoters <- max_states




##### ATP12A p-values from UNC lung association (published GWAS on the Strug website)
##### Will take and plot the fixed meta-GWAS p-values
##### Taking the column ALL rather than FF (FF are the deltaF508 homozygous-contrained individuals' association)

lung_assoc_atp12a <- read.table("atp12a_lung_association_summary_stats.txt", header=T, stringsAsFactors = F)
lung_assoc_atp12a <- lung_assoc_atp12a[,c(1,3,2,10)]
names(lung_assoc_atp12a) <- c("chr","snp","pos","p")
lung_assoc_atp12a <- cbind(lung_assoc_atp12a, MarkerName=paste0("chr", lung_assoc_atp12a$chr,":",lung_assoc_atp12a$pos))


#HNE_eqtl <- read.table("/Users/naim/Documents/Strug/CF/eqtl/Nasal_RNAseq/SLC6A14_HNE_eqtl_by_Gengming.txt", header=T, stringsAsFactors = F)
HNE_eqtl <- read.table("/Users/naim/Documents/Strug/CF/eqtl/Nasal_RNAseq/all_nominal_SLC6A14_hg19_Cheng.txt", header=T, stringsAsFactors = F)
slc6a14_start <- 115567790
HNE_eqtl <- cbind(chr="X", HNE_eqtl)
HNE_eqtl <- cbind(HNE_eqtl, pos=HNE_eqtl$distance+slc6a14_start)
#colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "POS")] <- "pos"
#colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "nominal_p")] <- "pval_nominal"
colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "pval")] <- "pval_nominal"
HNE_eqtl <- cbind(MarkerName=paste0(HNE_eqtl$chr,":",HNE_eqtl$pos), HNE_eqtl)
#colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "gene")] <- "gene_id"
colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "gene2")] <- "gene_id"
HNE_eqtl <- HNE_eqtl[order(HNE_eqtl$pos),]
write.table(HNE_eqtl, paste0(gtex_dir, "HNE_slc6a14_eQTLs_mod.txt"), quote=F, col.names=T, row.names=F)


# HNE males and females 6a14 eQTL data from Cheng:
HNE_eqtl_males <- read.table("/Users/naim/Documents/Strug/CF/eqtl/Nasal_RNAseq/male_nominal_SLC6A14_hg19_Cheng.txt", header=T, stringsAsFactors = F)
HNE_eqtl_females <- read.table("/Users/naim/Documents/Strug/CF/eqtl/Nasal_RNAseq/female_nominal_SLC6A14_hg19_Cheng.txt", header=T, stringsAsFactors = F)

fix_HNE_file <- function(HNE_eqtl) {
  slc6a14_start <- 115567790
  HNE_eqtl <- cbind(chr="X", HNE_eqtl)
  HNE_eqtl <- cbind(HNE_eqtl, pos=HNE_eqtl$distance+slc6a14_start)
  colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "pval")] <- "pval_nominal"
  HNE_eqtl <- cbind(MarkerName=paste0(HNE_eqtl$chr,":",HNE_eqtl$pos), HNE_eqtl)
  colnames(HNE_eqtl)[which(colnames(HNE_eqtl) %in% "gene2")] <- "gene_id"
  HNE_eqtl <- HNE_eqtl[order(HNE_eqtl$pos),]
}

HNE_eqtl_males <- fix_HNE_file(HNE_eqtl_males)
HNE_eqtl_females <- fix_HNE_file(HNE_eqtl_females)
write.table(HNE_eqtl_males, paste0(gtex_dir, "HNE_slc6a14_eQTLs_males_mod.txt"), quote=F, col.names=T, row.names=F)
write.table(HNE_eqtl_females, paste0(gtex_dir, "HNE_slc6a14_eQTLs_females_mod.txt"), quote=F, col.names=T, row.names=F)



#========================================= Functions ======================================================
plot_correlations <- function(genename, chr, tissue) {
  gene <- tolower(genename)
  eqtls <- read.table(paste0(gtex_dir,tissue, "_",gene,"_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
  gene_tissue_merge <- merge(eqtls, MI_suggestive_pvalues,by="MarkerName")
  if (length(names(gene_tissue_merge)[which(names(gene_tissue_merge) %in% "p.value")]) != 0) {
    names(gene_tissue_merge)[which(names(gene_tissue_merge) %in% "p.value")] <- "pval_nominal"
  }
  
  xrange <- c(min(gene_tissue_merge$V3), max(gene_tissue_merge$V3))
  yrange <- c(0, max(c(-log10(gene_tissue_merge$V4), -log10(gene_tissue_merge$pval_nominal)))+3)
  
  pdf(paste0(toupper(genename),"_GWAS_vs_GTEx_",tissue,"_eQTL_pvalues.pdf"))
  plot(xrange,yrange,type="n", main=paste(toupper(genename),"GWAS and GTEx",tissue, "eQTL P-values"), 
       xlab=paste0("chr",chr," coordinates"), ylab="-log10P-values")
  lines(gene_tissue_merge$V3,-log10(gene_tissue_merge$V4),col="red")
  points(gene_tissue_merge$V3,-log10(gene_tissue_merge$V4), col="red")
  lines(x=gene_tissue_merge$V3,-log10(gene_tissue_merge$pval_nominal))
  points(x=gene_tissue_merge$V3,-log10(gene_tissue_merge$pval_nominal))
  legend("topleft",c("MI GWAS p-values", paste("GTEx",tissue,"eQTL p-values")), col=c("red","black"), lty=1)
  dev.off()
  
  # P-P plot
  pdf(paste0(toupper(genename),"_",tissue,"_eQTL_P_vs_GWAS_P.pdf"))
  plot(-log10(gene_tissue_merge$V4), -log10(gene_tissue_merge$pval_nominal), 
       main=paste(toupper(genename),tissue,"P-P plot"),xlab="GWAS -log10P",ylab=paste(tissue,"eQTL -log10P"))
  rho=cor(-log10(gene_tissue_merge$V4), -log10(gene_tissue_merge$pval_nominal))
  y=-log10(gene_tissue_merge$pval_nominal)
  x=-log10(gene_tissue_merge$V4)
  abline(a = coef(lm(y ~ x))[1], b=coef(lm(y ~ x))[2],col="red")
  mtext(paste0("correlation = ", signif(rho,4)), side=3)
  dev.off()
}
#-----------------------------------------------------------------
smoothing <- function(x,y,xrange,yrange,window_size) {
  window_partition <- window_size
  windowing <- (xrange[2]-xrange[1])/window_partition
  curr <- xrange[1]
  smooth_curve_x <- NULL
  smooth_curve_y <- NULL
  while(curr < min(xrange[2],x[length(x)]) & length(which(x>=curr & x<=(curr+windowing)))>0) {
    indices <- which(x>=curr & x<=(curr+windowing))
    smooth_curve_x <- c(smooth_curve_x, mean(x[indices]))
    smooth_curve_y <- c(smooth_curve_y, max(y[indices]))
    curr <- curr + windowing + 1
  }
  return(list(smooth_curve_x, smooth_curve_y))
}
#-----------------------------------------------------------------
smoothing2 <- function(x,y,xrange,yrange,window_size) {
  window_partition <- window_size
  windowing <- (xrange[2]-xrange[1])/window_partition
  curr <- xrange[1]
  smooth_curve_x <- NULL
  smooth_curve_y <- NULL
  while(curr < min(xrange[2],x[length(x)]) & length(which(x>=curr & x<=(curr+windowing)))>0) {
    indices <- which(x>=curr & x<=(curr+windowing))
    smooth_curve_x <- c(smooth_curve_x, subset(x[indices], y[indices] %in% max(y[indices]))[1])
    smooth_curve_y <- c(smooth_curve_y, max(y[indices]))
    curr <- curr + windowing + 1
  }
  return(list(smooth_curve_x, smooth_curve_y))
}
#-----------------------------------------------------------------
plot_eQTL <- function(genename, tissues, chr, start, end, ymax, window_size, 
                      draw_dots=FALSE, draw_linefit=TRUE,
                      title="", lty=1, colour="black", superimpose=FALSE, yaxis="left") {
  tissue_count <- 0
  yaxis_drawn <- FALSE
  something_drawn_already <- FALSE
  for(i in tissues) {
    tissue_count <- tissue_count + 1
    gene_tissue_eqtls <- read.table(paste0(gtex_dir,i,"_",tolower(genename), "_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
    if (length(names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")]) != 0) {
      names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")] <- "pval_nominal"
    }
    if(nrow(gene_tissue_eqtls)!=0) {
      
      x <- as.numeric(gsub(paste0(chr,":(\\d*?)"),"\\1",gene_tissue_eqtls[,1]))
      y <- -log10(gene_tissue_eqtls$pval_nominal)
      
      xrange <- c(start,end)
      #yrange <- c(0,max(y[min(which(x>=start)):max(which(x<=end))]))
      yrange <- c(-2,ymax+0.5)
      smooth_x_y <- smoothing2(x,y,xrange,yrange,window_size)
      
      # smooth_x_y <- smooth.spline(x,y,df=2000)
      # lo <- loess(y~x, degree=2)
      # xl <- seq(min(x),max(x),(max(x)-min(x))/10000)
      # lines(xl/10^6, predict(lo,xl),col="green",lwd=2)
      
      if(!something_drawn_already & !superimpose) {
        plot(xrange/10^6, yrange, main=title, 
             xlab=paste0("chr",chr, " coordinates (Mb)"), ylab="GTEx eQTL -log10(P-values)", type="n", bty="L")
        something_drawn_already <- TRUE
      } else {
        par(new=T)
        plot(xrange/10^6, yrange,
             type="n", bty="L",
             xaxt='n', yaxt='n', ann=FALSE)
        if(yaxis %in% "right" & !yaxis_drawn) {
          par(xpd=T)
          par(las=3)
          axis(side=4)
          mtext("GTEx eQTL -log10(P-values)", side=4, line=2)
          yaxis_drawn <- TRUE
        }
        if(draw_dots) points(x/10^6,y,col="black",pch=16,cex=0.75, col=colour[tissue_count])
        if(draw_linefit) lines(smooth_x_y[[1]]/10^6, smooth_x_y[[2]], lwd=2, lty=lty[tissue_count], col=colour[tissue_count])
      }
    }
  }
}
#-----------------------------------------------------------------
plot_eQTL2 <- function(eQTL_filenames, ENSGname, chr, start, end, ymax, window_size, 
                      draw_dots=FALSE, draw_linefit=TRUE,
                      title="", lty=1, colour="black", superimpose=FALSE, yaxis="left") {
  # eQTL_filename columns example:
  #MarkerName           gene_id           variant_id tss_distance pval_nominal      slope slope_se
  #X:114569348 ENSG00000087916.7  X_114569348_G_A_b37      -998442    0.3382370  0.1395100 0.145056
  #X:114569770 ENSG00000087916.7  X_114569770_C_T_b37      -998020    0.3435240  0.1364990 0.143494
  
  # MarkerName must be the first column
  # ensure file is small enough to be read by read.table
  # pval_nominal column name is important to be named either pval_nominal or p.value
  # gene_id column is needed to subset the region for the p-values for eQTL of GOI
  # all other columns are not used
  
  tissue_count <- 0
  yaxis_drawn <- FALSE
  something_drawn_already <- FALSE
  for(i in eQTL_filenames) {
    tissue_count <- tissue_count + 1
    gene_tissue_eqtls <- read.table(i, header=T, stringsAsFactors = F)
    if (length(names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")]) != 0) {
      names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")] <- "pval_nominal"
    }
    chrom <- as.character(gsub("(.*):(.*)", "\\1", gene_tissue_eqtls$MarkerName))
    pos <- as.integer(gsub("(.*):(.*)", "\\2", gene_tissue_eqtls$MarkerName))
    gene_tissue_eqtls <- subset(gene_tissue_eqtls, gene_tissue_eqtls$gene_id %in% as.character(ENSGname) &
                                  chrom %in% as.character(chr) & pos >= start & pos <= end)
    if(nrow(gene_tissue_eqtls)!=0) {
      
      x <- as.numeric(gsub(paste0(chr,":(\\d*?)"),"\\1",gene_tissue_eqtls[,1]))
      y <- -log10(gene_tissue_eqtls$pval_nominal)
      
      xrange <- c(start,end)
      #yrange <- c(0,max(y[min(which(x>=start)):max(which(x<=end))]))
      yrange <- c(-2,ymax+0.5)
      smooth_x_y <- smoothing2(x,y,xrange,yrange,window_size)
      
      # smooth_x_y <- smooth.spline(x,y,df=2000)
      # lo <- loess(y~x, degree=2)
      # xl <- seq(min(x),max(x),(max(x)-min(x))/10000)
      # lines(xl/10^6, predict(lo,xl),col="green",lwd=2)
      
      if(!something_drawn_already & !superimpose) {
        plot(xrange/10^6, yrange, main=title, 
             xlab=paste0("chr",chr, " coordinates (Mb)"), ylab="GTEx eQTL -log10(P-values)", type="n", bty="L")
        something_drawn_already <- TRUE
      } else {
        par(new=T)
        plot(xrange/10^6, yrange,
             type="n", bty="L",
             xaxt='n', yaxt='n', ann=FALSE)
        if(yaxis %in% "right" & !yaxis_drawn) {
          par(xpd=T)
          par(las=3)
          axis(side=4)
          mtext("GTEx eQTL -log10(P-values)", side=4, line=2)
          yaxis_drawn <- TRUE
        }
        if(draw_dots) points(x/10^6,y,col="black",pch=16,cex=0.75, col=colour[tissue_count])
        if(draw_linefit) lines(smooth_x_y[[1]]/10^6, smooth_x_y[[2]], lwd=2, lty=lty[tissue_count], col=colour[tissue_count])
      }
    }
  }
}
#-----------------------------------------------------------------
get_eQTL_yrange <- function(genename, tissues) {
  ymax <- 0
  for(i in tissues) {
    gene_tissue_eqtls <- read.table(paste0(gtex_dir,i,"_",tolower(genename), "_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
    if (length(names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")]) != 0) {
      names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")] <- "pval_nominal"
    }
    if(nrow(gene_tissue_eqtls)!=0) {
      ymax2 <- max(-log10(gene_tissue_eqtls$pval_nominal))
      ymax <- ifelse(ymax2>ymax,ymax2,ymax) 
    }
  }
  return(ymax)
}
#-----------------------------------------------------------------
get_eQTL_yrange2 <- function(eQTL_filenames, ENSGname, chr, startbp, endbp) {
  ymax <- 0
  for(i in eQTL_filenames) {
    gene_tissue_eqtls <- read.table(i, header=T, stringsAsFactors = F)
    chrom <- as.character(gsub("(.*):(.*)", "\\1", gene_tissue_eqtls$MarkerName))
    pos <- as.integer(gsub("(.*):(.*)", "\\2", gene_tissue_eqtls$MarkerName))
    gene_tissue_eqtls <- subset(gene_tissue_eqtls, gene_tissue_eqtls$gene_id %in% as.character(ENSGname) & 
                                  chrom %in% as.character(chr) & pos >= startbp & pos <= endbp)
    if (length(names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")]) != 0) {
      names(gene_tissue_eqtls)[which(names(gene_tissue_eqtls) %in% "p.value")] <- "pval_nominal"
    }
    if(nrow(gene_tissue_eqtls)!=0) {
      ymax2 <- max(-log10(gene_tissue_eqtls$pval_nominal))
      ymax <- ifelse(ymax2>ymax,ymax2,ymax) 
    }
  }
  return(ymax)
}
#-----------------------------------------------------------------
plot_locuszoom <- function(gwasResult, chrom, start, end, yrange, ld_filename, lead_snpname, transparency, superimpose=TRUE, yaxis="left",
                            r2colors=c("purple","red","orange","green","cyan","blue","black")) {
   if(chrom=="X") chrom=23
   if(gwasResult$chr[1]=="X") {
     gwasResult$chr <- gsub("X","23", gwasResult$chr)
     gwasResult$MarkerName <- gsub("(chrX:|X:)\\d*?","chr23:",gwasResult$MarkerName)
   }
   GWAS_pvalues <- subset(gwasResult, (gwasResult$chr %in% chrom) &
                            (gwasResult$pos >= start) & (gwasResult$pos <= end))
   if(length(which(GWAS_pvalues$p==0))!=0) { # if there are any p-values with a 0 value
    GWAS_pvalues[which(GWAS_pvalues$p==0),]$p <- 2.2e-16  # set them to the lowest p-value possible by R
  }
  ld <- read.table(ld_filename,
                   header=T, stringsAsFactors = F)
  names(ld)[which(names(ld) %in% "SNP_B")] <- "MarkerName"
  ld_mergedset <- merge(GWAS_pvalues, ld[,c(6,7)], by="MarkerName")
  
  lead_snp <- subset(ld_mergedset, ld_mergedset$MarkerName %in% lead_snpname)
  r2_80_group <- subset(ld_mergedset, ld_mergedset$R2>=0.8)
  r2_60_group <- subset(ld_mergedset, ld_mergedset$R2>=0.6 & ld_mergedset$R2<0.8)
  r2_40_group <- subset(ld_mergedset, ld_mergedset$R2>=0.4 & ld_mergedset$R2<0.6)
  r2_20_group <- subset(ld_mergedset, ld_mergedset$R2>=0.2 & ld_mergedset$R2<0.4)
  r2_lt_20_group <- subset(ld_mergedset, ld_mergedset$R2<0.2)
  no_ld_info_snps <- subset(GWAS_pvalues, !(GWAS_pvalues$MarkerName %in% ld_mergedset$MarkerName))
  
  lowest_p <- min(ld_mergedset$p)
  
  if(superimpose==T & yaxis %in% "right") {
    plot(c(start/10^6, end/10^6), c(yrange[1],max(c(-log10(lowest_p),yrange[2]))), axes=F, xlab=NA, ylab=NA, type="n")
    axis(side=4)
    mtext("GWAS -log10(P-values)", side=4, line=2)
    points(no_ld_info_snps$pos/10^6, -log10(no_ld_info_snps$p), pch=21, bg=alpha(r2colors[7],transparency),cex=1.5,lwd=0.25)
    points(r2_lt_20_group$pos/10^6, -log10(r2_lt_20_group$p), pch=21, bg=alpha(r2colors[6],transparency),cex=1.5,lwd=0.25)
    points(r2_20_group$pos/10^6, -log10(r2_20_group$p), pch=21, bg=alpha(r2colors[5],transparency),cex=1.5,lwd=0.25)
    points(r2_40_group$pos/10^6, -log10(r2_40_group$p), pch=21, bg=alpha(r2colors[4],transparency),cex=1.5,lwd=0.25)
    points(r2_60_group$pos/10^6, -log10(r2_60_group$p), pch=21, bg=alpha(r2colors[3],transparency),cex=1.5,lwd=0.25)
    points(r2_80_group$pos/10^6, -log10(r2_80_group$p), pch=21, bg=alpha(r2colors[2],transparency),cex=1.5,lwd=0.25)
    points(lead_snp$pos/10^6, -log10(lead_snp$p), pch=23, col=alpha(r2colors[1],transparency), bg=alpha("purple",transparency),cex=1.5,lwd=0.25)
    text(x=lead_snp$pos/10^6,y=-log10(lead_snp$p)+1.5, labels=lead_snp$snp)
  } else if (superimpose==T & yaxis %in% "left") {
    points(no_ld_info_snps$pos/10^6, -log10(no_ld_info_snps$p), pch=21, bg=alpha(r2colors[7],transparency),cex=1.5,lwd=0.25)
    points(r2_lt_20_group$pos/10^6, -log10(r2_lt_20_group$p), pch=21, bg=alpha(r2colors[6],transparency),cex=1.5,lwd=0.25)
    points(r2_20_group$pos/10^6, -log10(r2_20_group$p), pch=21, bg=alpha(r2colors[5],transparency),cex=1.5,lwd=0.25)
    points(r2_40_group$pos/10^6, -log10(r2_40_group$p), pch=21, bg=alpha(r2colors[4],transparency),cex=1.5,lwd=0.25)
    points(r2_60_group$pos/10^6, -log10(r2_60_group$p), pch=21, bg=alpha(r2colors[3],transparency),cex=1.5,lwd=0.25)
    points(r2_80_group$pos/10^6, -log10(r2_80_group$p), pch=21, bg=alpha(r2colors[2],transparency),cex=1.5,lwd=0.25)
    points(lead_snp$pos/10^6, -log10(lead_snp$p), pch=23, col=alpha(r2colors[1],transparency), bg=alpha("purple",transparency),cex=1.5,lwd=0.25)
    text(x=lead_snp$pos/10^6,y=-log10(lead_snp$p)+1.5, labels=lead_snp$snp)
  } else {
    plot(x=c(start/10^6, end/10^6), y=c(yrange[1],max(c(-log10(lowest_p),yrange[2]))), 
         xlab=paste0("chr",chrom," coordinates (Mb)"), ylab="GWAS -log10(P-values)",type="n")
    if(yaxis %in% "right") {
      axis(side=4)
      mtext("GWAS -log10(P-values)", side=4, line=2)
    }
    points(no_ld_info_snps$pos/10^6, -log10(no_ld_info_snps$p), pch=21, bg=alpha(r2colors[7],transparency),cex=1.5,lwd=0.25)
    points(r2_lt_20_group$pos/10^6, -log10(r2_lt_20_group$p), pch=21, bg=alpha(r2colors[6],transparency),cex=1.5,lwd=0.25)
    points(r2_20_group$pos/10^6, -log10(r2_20_group$p), pch=21, bg=alpha(r2colors[5],transparency),cex=1.5,lwd=0.25)
    points(r2_40_group$pos/10^6, -log10(r2_40_group$p), pch=21, bg=alpha(r2colors[4],transparency),cex=1.5,lwd=0.25)
    points(r2_60_group$pos/10^6, -log10(r2_60_group$p), pch=21, bg=alpha(r2colors[3],transparency),cex=1.5,lwd=0.25)
    points(r2_80_group$pos/10^6, -log10(r2_80_group$p), pch=21, bg=alpha(r2colors[2],transparency),cex=1.5,lwd=0.25)
    points(lead_snp$pos/10^6, -log10(lead_snp$p), pch=23, col=alpha(r2colors[1],transparency), bg=alpha("purple",transparency),cex=1.5,lwd=0.25)
    text(x=lead_snp$pos/10^6,y=-log10(lead_snp$p)+1.5, labels=lead_snp$snp)
  }
}
#-----------------------------------------------------------------
plot_r2_legend <- function(leadSnpColor="purple",r2color80="red", r2color60="orange", r2color40="green", 
                           r2color20="cyan", r2color0="blue", r2colorNA="black") {
  num_categories <- 7
  topy <- 2
  bottomy <- 0
  margin <- 0.1
  textsize <- 1
  circlesize <- 2
  
  plot(c(0,0.5),c(bottomy-margin,topy+margin),type='n',xlab='',ylab='', axes = F)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*0,pch=21, bg=r2colorNA,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*0,expression('N/A'),pos=4, pch=21,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*1,pch=21, bg=r2color0,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*1,expression('<0.2'),pos=4, pch=21,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*2,pch=21, bg=r2color20,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*2,expression('0.2'),pos=4,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*3,pch=21, bg=r2color40,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*3,expression('0.4'),pos=4,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*4,pch=21, bg=r2color60,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*4,expression('0.6'),pos=4,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*5,pch=21, bg=r2color80,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*5,expression('0.8'),pos=4,cex=textsize)
  
  points(0.1,bottomy + ((topy-bottomy)/(num_categories-1))*6,pch=23, bg=leadSnpColor,cex=circlesize+1.5)
  text(0.11,bottomy + ((topy-bottomy)/(num_categories-1))*6,expression('REF SNP'),pos=4, pch=21,cex=textsize)
}
#-----------------------------------------------------------------
draw_gene <- function(y, exon_starts, exon_ends, genename, window_min, window_max, ...) {
  exon_height <- 0.5
  intron_height <- 0.05
  
  trunc_start <- which(exon_starts > window_min)[1] != 1
  trunc_end <- exon_ends[length(exon_ends)] > window_max
  new_exon_starts <- exon_starts
  new_exon_ends <- exon_ends
  if(trunc_start & trunc_end) {
    new_exon_starts <- c(exon_starts[which(exon_starts >= window_min & exon_starts <= window_max)], window_max)
    new_exon_ends <- c(window_min, exon_ends[which(exon_ends >= window_min & exon_ends <= window_max)])
    gene_length <- length(new_exon_starts)
    rect(new_exon_starts[1:(length(new_exon_starts)-1)], y, new_exon_ends[2:length(new_exon_ends)], y+exon_height, col="black")
    rect(new_exon_ends[1:gene_length], y+(exon_height/2)-intron_height, 
         new_exon_starts[1:gene_length], y+(exon_height/2)+intron_height,col="black")
  }
  else if(trunc_start) {
    new_exon_starts <- c(exon_starts[which(exon_starts >= window_min & exon_starts <= window_max)])
    new_exon_ends <- c(window_min, exon_ends[which(exon_ends >= window_min & exon_ends <= window_max)])
    gene_length <- length(new_exon_starts)
    rect(new_exon_starts[1:(length(new_exon_starts))], y, new_exon_ends[2:length(new_exon_ends)], y+exon_height, col="black")
    rect(new_exon_ends[1:gene_length], y+(exon_height/2)-intron_height, 
         new_exon_starts[1:gene_length], y+(exon_height/2)+intron_height,col="black")
  }
  else if(trunc_end) {
    new_exon_starts <- c(exon_starts[which(exon_starts >= window_min & exon_starts <= window_max)], window_max)
    new_exon_ends <- c(exon_ends[which(exon_ends >= window_min & exon_ends <= window_max)])
    gene_length <- length(new_exon_starts)
    rect(new_exon_starts[1:(length(new_exon_starts)-1)], y, new_exon_ends[1:length(new_exon_ends)], y+exon_height, col="black")
    rect(new_exon_ends[1:length(new_exon_ends)], y+(exon_height/2)-intron_height, 
         new_exon_starts[2:length(new_exon_starts)], y+(exon_height/2)+intron_height,col="black")
  } else {
    gene_length <- length(new_exon_starts)
    rect(exon_starts, y, exon_ends, y+exon_height, col="black")
    rect(exon_ends[1:(gene_length-1)], y+(exon_height/2)-intron_height,
         exon_starts[2:gene_length], y+(exon_height/2)+intron_height,col="black")
  }
  # exon_starts <- new_exon_starts
  # exon_ends <- new_exon_ends
  text(new_exon_ends[length(new_exon_ends)],y+exon_height/2,genename,col="red",pos=4, ...)
}
#-----------------------------------------------------------------
colour_regions_track <- function(chr, coord_range, ystart, trackheight, positive_regions, positive_regions_col, units="Mb") {
  # PRE: chr is the chromosome integer that will be plot; for chrX, enter 23
  #      coord_range is the x-axis coordinate range of interest that we want to display; eg. c(205860000, 205940000)
  #      newplot is boolean to whether a new plot should be drawn (specify TRUE), or add to existing plot (specify FALSE)
  #      ystart is the y coordinate of the top of the track
  #      trackheight is the desired height of the track
  #      positive_regions specifies the regions that have positive hits (eg. enhancer regions); format as chr1:205889001-205889026
  #      positive_regions_col specifies the fill colour of the regions with positive hits
  #      units can be Mb or kb or bp
  # POST: a track is drawn with the positive regions coloured
  
  positive_regions <- cbind(positive_regions, CHR=gsub("chr(\\d*|X):\\d*-\\d*","\\1",positive_regions))
  positive_regions <- cbind(positive_regions, START=gsub("chr(\\d*|X):(\\d*)-(\\d*)","\\2",positive_regions[,1]))
  positive_regions <- cbind(positive_regions, END=gsub("chr(\\d*|X):(\\d*)-(\\d*)","\\3",positive_regions[,1]))
  positive_regions[,'CHR'] <- ifelse(positive_regions[,'CHR'] %in% "X", 23, positive_regions[,'CHR'])
  positive_regions2 <- as.matrix(cbind(CHR=as.numeric(positive_regions[,2]), START=as.numeric(positive_regions[,3]), END=as.numeric(positive_regions[,4])))
  
  positive_regions3 <- subset(positive_regions2, positive_regions2[,'CHR'] %in% chr & positive_regions2[,'START'] >= coord_range[1] & positive_regions2[,'END'] <= coord_range[2])
  divisor <- 1
  if(units %in% c("Mb","MB","mb")) {
    divisor <- 10^6
  } else if(units %in% c("Kb","KB","kb")) {
    divisor <- 10^3
  }
  
  if(nrow(positive_regions3) != 0) {
    rect(positive_regions3[, 'START']/divisor, ystart - trackheight, positive_regions3[, 'END']/divisor, ystart, col=positive_regions_col, lwd=0.0001)
  }
  
}
#-----------------------------------------------------------------



#========================================= Main ======================================================

##==== eQTL Plots =======
transparency <- 1.0
transparency2 <- 1.0
tissues <- c("Pancreas", "Small_intestine", "Stomach", "Colon_Transverse", "Colon_Sigmoid")
cols <- c("black","blue","orange","green","red")
# tissues <- c("Pancreas", "Lung", "HNE")
# cols <- c("black","blue","orange")


#######################################################################################################
###########################################  Fig. 2A - SLC6A14  #######################################
#######################################################################################################

### SLC6A14 overall ####

chr <- "X"
bpstart <- 115200000
bpend <- 115700000

slc6a14_eQTL_files <- c("fastQTL_runs/Pancreas_FastQTL_X_114500000_116600000_no_age_fixed_slopeSE_eQTLs_mod.txt.gz"
                        ,"Small_Intestine_Terminal_Ileum.v7.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Stomach.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Colon_Transverse.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Colon_Sigmoid.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
)

pdf("v7_plots/SLC6A14_vs_eQTLs_overall_GTExV7.pdf")
#pdf("SLC6A14_locuszoom2.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+0,2.1+2)) # par()$mar
#par(mar=c(5.1,4.1,15.1,2.1)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs3788766_slc6a14.ld"
lead_snpname <- "chr23:115566839"
extra_y_height <- 5

ymax <- max(-log10(min(All_GWAS_pvalues$p))) + extra_y_height
#ymax=16
#ymin <- (-1)*(2/(get_eQTL_yrange("SLC6A14",tissues)+1))*ymax
ymin <- -4
plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
               r2colors=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"))
# plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
#                r2colors=c("purple","red","orange","green","cyan","blue","black"))

#plot_r2_legend("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "black")
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7879546_agtr2.ld"
lead_snpname <- as.character(lung_assoc_overall_UNC$MarkerName[which(lung_assoc_overall_UNC$p == min(lung_assoc_overall_UNC$p))])
lead_snpname <- gsub("chrX", "chr23", lead_snpname)
plot_locuszoom(lung_assoc_overall_UNC, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = T, yaxis="left",
               r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))
#plot_r2_legend("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","black")
# plot_eQTL(genename="SLC6A14", tissues,chr="X",start=bpstart,end=bpend,window_size=25,ymax=get_eQTL_yrange("SLC6A14",tissues)+1, superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="SLC6A14",lty=seq(1:length(tissues)), colour=cols)
# legend("top", tissues, title="SLC6A14 eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
plot_eQTL2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7",chr="X",start=bpstart,end=bpend,window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7", chr, bpstart, bpend)+1, 
           superimpose=T, yaxis="right",
          draw_dots = F, draw_linefit = T, title="SLC6A14",lty=seq(1:length(tissues)), colour=cols)
legend("top", tissues, title="SLC6A14 eQTL P-values",
       lty=rep(1,3), lwd=rep(2,length(tissues)),col=cols, bty='n')


legend("topleft", title=expression("Lung GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"), bty='n')
legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

# genename="slc6a14"
# gene_tissue_eqtls <- read.table(paste0(gtex_dir,i,"_",tolower(genename), "_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
# ymax=get_eQTL_yrange("SLC6A14",tissues)+1
# pancreas_eqtl <- data.frame(snp=as.character(gene_tissue_eqtls$MarkerName), chr="X", pos=as.integer(gsub("X:(\\d*?)","\\1",gene_tissue_eqtls$MarkerName)),
#                             p=gene_tissue_eqtls$pval_nominal, MarkerName=gene_tissue_eqtls$MarkerName, stringsAsFactors = F)
# plot_locuszoom(pancreas_eqtl, "X", 115400000,115700000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = F, yaxis="left",
#                r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))


dumbo = 0

agtr2_exon_starts <- c(115301957,115302220,115303498)/10^6
agtr2_exon_ends <- c(115302069,115302280,115306225)/10^6
draw_gene(y=-1.5+dumbo, agtr2_exon_starts,agtr2_exon_ends, "AGTR2 -->",bpstart/10^6,bpend/10^6)

slc6a14_exon_starts <- c(115567746,115568957,115572133,115573854,115574810,115576085,115577906,115582606,115584181,115585489,115586142,115586522,115588774,115589974)/10^6
slc6a14_exon_ends <- c(115567925,115569123,115572265,115574016,115574958,115576218,115578047,115582835,115584307,115585608,115586242,115586632,115588942,115592625)/10^6
draw_gene(y=-1.25+dumbo,slc6a14_exon_starts,slc6a14_exon_ends,"SLC6A14 -->",bpstart/10^6,bpend/10^6)

ct83_exon_starts <- c(115592852,115593942)/10^6
ct83_exon_ends <- c(115593174,115594194)/10^6
ct83_length <- length(ct83_exon_starts)
draw_gene(y=-1.85+dumbo,ct83_exon_starts,ct83_exon_ends,"<-- CT83",bpstart/10^6,bpend/10^6)

#dev.off()



####### Where are the enhancers and promoters? ########
### None to be found in any of the tissues for chrX


ystart = 6.5
step = 0.5
trackheight = 0.4
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=23, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight,
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=23, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight,
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(bpend/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()



#######################################################################################################
###########################################  Fig. 2B - SLC26A9  #######################################
#######################################################################################################

chr <- "1"
bpstart <- 205860000
bpend <- 205940000

slc26a9_eQTL_files <- c("Pancreas.allpairs_fixed_chr1_205000000_208000000_eQTLs_mod.txt.gz"
                        ,"Small_Intestine_Terminal_Ileum.v7.allpairs_fixed_chr1_205000000_208000000_eQTLs_mod.txt.gz"
                        ,"Stomach.allpairs_fixed_chr1_205000000_208000000_eQTLs_mod.txt.gz"
                        ,"Colon_Transverse.allpairs_fixed_chr1_205000000_208000000_eQTLs_mod.txt.gz"
                        ,"Colon_Sigmoid.allpairs_fixed_chr1_205000000_208000000_eQTLs_mod.txt.gz"
)

pdf("v7_plots/SLC26A9_vs_eQTLs_GTExV7.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+0,2.1+2)) # par()$mar
#par(mar=c(5.1,4.1,4.1,2.1)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7549173_slc26a9.ld"
lead_snpname <- "chr1:205906897"
extra_y_height <- 7
region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 1 & All_GWAS_pvalues$pos >= 205860000 & All_GWAS_pvalues$pos <= 205940000)
ymax <- max(-log10(min(region$p))) + extra_y_height
ymin <- (-1)*(2/(get_eQTL_yrange2(paste0(gtex_dir, slc26a9_eQTL_files), "ENSG00000174502.14", chr, bpstart, bpend)+1))*ymax
#ymin <- (-1)*(2/(10+1))*ymax
plot_locuszoom(All_GWAS_pvalues, "1",bpstart,bpend, yrange=c(ymin,ymax), ld_filename, lead_snpname, transparency, superimpose=F, yaxis="left")
plot_eQTL2(paste0(gtex_dir, slc26a9_eQTL_files), "ENSG00000174502.14","1",bpstart,bpend,window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir, slc26a9_eQTL_files), "ENSG00000174502.14", chr, bpstart, bpend),  
           superimpose=T, yaxis="right",
           draw_dots = F, draw_linefit = T, title="SLC26A9",lty=seq(1:length(tissues)), colour=cols)
legend("topleft", tissues, title="SLC26A9 Gene Expression eQTL P-values",
       lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

LOC284581_exon_starts <- c(205831206,205860283,205863295)/10^6
LOC284581_exon_ends <- c(205831726,205860437,205865215)/10^6
SLC26A9_exon_starts <- c(205882176,205886410,205887967,205889303,205890693,205892209,205892462,205892671,205893510,205895662,205896338,205896619,205897029,205897954,205898331,205899019,205900987,205901829,205902072,205904823,205912492)/10^6
SLC26A9_exon_ends <- c(205884532,205886482,205888113,205889358,205890975,205892323,205892555,205892741,205893617,205895758,205896416,205896733,205897177,205898037,205898484,205899184,205901163,205901940,205902212,205904966,205912588)/10^6

draw_gene(y=-1, LOC284581_exon_starts, LOC284581_exon_ends, "LOC284581 -->", (bpstart-2500)/10^6, bpend/10^6)
draw_gene(y=-1.75, SLC26A9_exon_starts, SLC26A9_exon_ends, "<-- SLC26A9", bpstart/10^6, bpend/10^6)

#dev.off()

# y=-1;exon_starts = SLC26A9_exon_starts; exon_ends = SLC26A9_exon_ends; genename="<-- SLC26A9"; window_min = bpstart/10^6; window_max = bpend/10^6
# y=-1;exon_starts = LOC284581_exon_starts; exon_ends = LOC284581_exon_ends; genename="LOC284581 -->"; window_min = bpstart/10^6; window_max = bpend/10^6

# ###### SLC41A1 eQTLs drawn instead ######
# plot_eQTL("SLC41A1", tissues, "1",205860000,205940000,window_size=25,ymax=get_eQTL_yrange("SLC26A9",tissues),  superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="SLC41A1",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="SLC41A1 Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
# 
# 
# ###### PM20D1 eQTLs drawn instead ######
# plot_eQTL("PM20D1", tissues, "1",205860000,205940000,window_size=25,ymax=10,  superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="PM20D1",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="PM20D1 Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
# 



####### Where are the enhancers and promoters? ########

ystart = 5.6
step = 0.25
trackheight = 0.2
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=1, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight, 
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=1, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight, 
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(205930000/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()




#######################################################################################################
#############################################  Fig 2C - ATP12A  #######################################
#######################################################################################################

chr <- "13"
bpstart <- 25200000
bpend <- 25350000

atp12a_eQTL_files <- c(
  "Pancreas.allpairs_fixed_chr13_24200000_26300000_eQTLs_mod.txt.gz"
  ,"Small_Intestine_Terminal_Ileum.v7.allpairs_fixed_chr13_24200000_26300000_eQTLs_mod.txt.gz"
  ,"Stomach.allpairs_fixed_chr13_24200000_26300000_eQTLs_mod.txt.gz"
  ,"fastQTL_runs/Colon_Transverse_FastQTL_13_24200000_26300000_no_age_fixed_slopeSE_eQTLs_mod.txt.gz"
  ,"Colon_Sigmoid.allpairs_fixed_chr13_24200000_26300000_eQTLs_mod.txt.gz"
)

pdf("v7_plots/ATP12A_vs_eQTLs_GTExV7.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+1,2.1+3)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs61948108_atp12a.ld"
lead_snpname <- "chr13:25282819"
extra_y_height <- 5
region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 13 & All_GWAS_pvalues$pos >= bpstart & All_GWAS_pvalues$pos <= bpend)
# ymax <- max(-log10(min(region$p))) + extra_y_height
ymax <- 15
#ymin <- (-1)*(2/(get_eQTL_yrange("ATP12A",tissues)+1))*ymax
ymin=-0.5
plot_locuszoom(All_GWAS_pvalues, "13",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")

# For overlay plot:
# plot_locuszoom(All_GWAS_pvalues, "13",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
#                r2colors=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"))

############################################################
# Code to overlay with lung association p-values from UNC (meta-GWAS fixed ALL (ie. Not just confined to deltaF508 homozygous))
############################################################
# ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs9553382_atp12a.ld"
# lead_snpname <- as.character(lung_assoc_atp12a$MarkerName[which(lung_assoc_atp12a$p[-56] == min(lung_assoc_atp12a$p[-56]))])  # lowest P is 13:25214757:T_TAT (no LD info for this one, so choose 2nd lowest P)
# plot_locuszoom(lung_assoc_atp12a[-56,], "13",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = T, yaxis="left",
#                r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))


############################################################
# Overlay condition ends
############################################################
plot_eQTL2(paste0(gtex_dir,atp12a_eQTL_files), ENSGname = "ENSG00000075673.7", chr = chr,start=bpstart, end=bpend, window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir,atp12a_eQTL_files), "ENSG00000075673.7", chr, bpstart, bpend), 
           superimpose=T, yaxis="right",
          draw_dots = F, draw_linefit = T, title="ATP12A",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="ATP12A Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')
legend("topleft", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"), bty='n')
# legend("topleft", title=expression("Lung GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"), bty='n')



# atp12a_pancreas_eqtls <- data.frame(snp=gene_tissue_eqtls$SNP, chr=13, pos=gsub("13:(\\d*?)","\\1",gene_tissue_eqtls$MarkerName), p=gene_tissue_eqtls$p.value, MarkerName=gene_tissue_eqtls$MarkerName, stringsAsFactors = F)
# which(atp12a_pancreas_eqtls$pos %in% 25282819)
# atp12a_pancreas_eqtls[4465,'snp'] <- "rs61948108"
# atp12a_pancreas_eqtls$MarkerName <- paste0("chr",atp12a_pancreas_eqtls$MarkerName)
# atp12a_pancreas_eqtls$pos <- as.integer(atp12a_pancreas_eqtls$pos)
# plot_locuszoom(atp12a_pancreas_eqtls, "13",25200000,25350000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")


TPTE2P6_exon_starts <- c(25154345,25157639,25160599,25161363,25163393,25168408,25169926,25171239,25171524)/10^6
TPTE2P6_exon_ends <- c(25155860,25157720,25160932,25161464,25163449,25168501,25169997,25171415,25171812)/10^6
ATP12A_exon_starts <- c(25254548,25255699,25259451,25262456,25263399,25264475,25264741,25265101,25266566,25266924,25268581,25272795,25274884,25276072,25280450,25281160,25281416,25283501,25283820,25284597,25284929,25285455,25285631)/10^6
ATP12A_exon_ends <- c(25254890,25255858,25259511,25262660,25263513,25264610,25264859,25265388,25266765,25267034,25268716,25272988,25275060,25276209,25280601,25281329,25281571,25283625,25283966,25284731,25285031,25285547,25285923)/10^6
RNF17_exon_starts <- c(25338300,25341409,25348950,25352432,25353804,25355981,25362125,25363485,25363835,25367179,25370274,25373532,25374503,25376518,25378425,25399756,25404621,25405995,25416178,25417881,25418789,25419098,25424478,25425593,25427992,25433138,25435405,25436850,25439010,25440281,25442737,25444708,25448251,25451134,25453324,25453874)/10^6
RNF17_exon_ends <- c(25338471,25341504,25349042,25352544,25353885,25356082,25362297,25363562,25363910,25367484,25370433,25373722,25374672,25376709,25378567,25399910,25404737,25406116,25416299,25418109,25418940,25419217,25424569,25425709,25428282,25433302,25435525,25436931,25439136,25440341,25442854,25444877,25448387,25451324,25453433,25454058)/10^6

dumbo = 32
draw_gene(y=-2.75+dumbo, TPTE2P6_exon_starts, TPTE2P6_exon_ends, "TPTE2P6 -->", bpstart/10^6, bpend/10^6)
draw_gene(y=-1.75+dumbo, ATP12A_exon_starts, ATP12A_exon_ends, "ATP12A -->", bpstart/10^6, bpend/10^6)
draw_gene(y=-3.5+dumbo, RNF17_exon_starts, RNF17_exon_ends, "RNF17 -->", bpstart/10^6, bpend/10^6)
#text(x=RNF17_exon_starts[1],y=-5+dumbo, "RNF17 -->",col="red")

# dev.off()



####### Where are the enhancers and promoters? ########

ystart = 41
step = 2.75
trackheight = 1
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=13, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight, 
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=13, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight, 
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(25360000/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()




#######################################################################################################
###########################################  Fig. 3A - SLC6A14  #######################################
#######################################################################################################

### SLC6A14 overall ####

chr <- "X"
bpstart <- 115200000
bpend <- 115700000

#cols[which(cols %in% "orange")] <- "darkgoldenrod3"

slc6a14_eQTL_files <- c("fastQTL_runs/Pancreas_FastQTL_X_114500000_116600000_no_age_fixed_slopeSE_eQTLs_mod.txt.gz"
                        ,"Lung.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"HNE_slc6a14_eQTLs_mod.txt"
)
tissues <- c("Pancreas", "Lung", "HNE")

pdf("v7_plots/SLC6A14_vs_eQTLs_overall_GTExV7_Lung_HNEs_Cheng_76HNEs.pdf")
#pdf("SLC6A14_locuszoom2.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+0,2.1+2)) # par()$mar
#par(mar=c(5.1,4.1,15.1,2.1)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs3788766_slc6a14.ld"
lead_snpname <- "chr23:115566839"
extra_y_height <- 5

ymax <- max(-log10(min(All_GWAS_pvalues$p))) + extra_y_height
#ymax=16
#ymin <- (-1)*(2/(get_eQTL_yrange("SLC6A14",tissues)+1))*ymax
ymin <- -4
plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
               r2colors=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"))
# plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
#                r2colors=c("purple","red","orange","green","cyan","blue","black"))

#plot_r2_legend("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "black")
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7879546_agtr2.ld"
lead_snpname <- as.character(lung_assoc_overall_UNC$MarkerName[which(lung_assoc_overall_UNC$p == min(lung_assoc_overall_UNC$p))])
lead_snpname <- gsub("chrX", "chr23", lead_snpname)
plot_locuszoom(lung_assoc_overall_UNC, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = T, yaxis="left",
               r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))
#plot_r2_legend("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","black")
# plot_eQTL(genename="SLC6A14", tissues,chr="X",start=bpstart,end=bpend,window_size=25,ymax=get_eQTL_yrange("SLC6A14",tissues)+1, superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="SLC6A14",lty=seq(1:length(tissues)), colour=cols)
# legend("top", tissues, title="SLC6A14 eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
plot_eQTL2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7",chr="X",start=bpstart,end=bpend,window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7", chr, bpstart, bpend)+1, 
           superimpose=T, yaxis="right",
           draw_dots = F, draw_linefit = T, title="SLC6A14",lty=c(1,2,4), colour=cols)
legend("top", tissues, title="SLC6A14 eQTL P-values",
       lty=c(1,2,4), lwd=rep(2,length(tissues)),col=cols, bty='n')


legend("topleft", title=expression("Lung GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"), bty='n')
legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
       pch=rep(16,7), col=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

# genename="slc6a14"
# gene_tissue_eqtls <- read.table(paste0(gtex_dir,i,"_",tolower(genename), "_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
# ymax=get_eQTL_yrange("SLC6A14",tissues)+1
# pancreas_eqtl <- data.frame(snp=as.character(gene_tissue_eqtls$MarkerName), chr="X", pos=as.integer(gsub("X:(\\d*?)","\\1",gene_tissue_eqtls$MarkerName)),
#                             p=gene_tissue_eqtls$pval_nominal, MarkerName=gene_tissue_eqtls$MarkerName, stringsAsFactors = F)
# plot_locuszoom(pancreas_eqtl, "X", 115400000,115700000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = F, yaxis="left",
#                r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))


dumbo = 0

agtr2_exon_starts <- c(115301957,115302220,115303498)/10^6
agtr2_exon_ends <- c(115302069,115302280,115306225)/10^6
draw_gene(y=-1.5+dumbo, agtr2_exon_starts,agtr2_exon_ends, "AGTR2 -->",bpstart/10^6,bpend/10^6)

slc6a14_exon_starts <- c(115567746,115568957,115572133,115573854,115574810,115576085,115577906,115582606,115584181,115585489,115586142,115586522,115588774,115589974)/10^6
slc6a14_exon_ends <- c(115567925,115569123,115572265,115574016,115574958,115576218,115578047,115582835,115584307,115585608,115586242,115586632,115588942,115592625)/10^6
draw_gene(y=-1.25+dumbo,slc6a14_exon_starts,slc6a14_exon_ends,"SLC6A14 -->",bpstart/10^6,bpend/10^6)

ct83_exon_starts <- c(115592852,115593942)/10^6
ct83_exon_ends <- c(115593174,115594194)/10^6
ct83_length <- length(ct83_exon_starts)
draw_gene(y=-1.85+dumbo,ct83_exon_starts,ct83_exon_ends,"<-- CT83",bpstart/10^6,bpend/10^6)

dev.off()









#######################################################################################################
############################################  AGTR2 - Fig. S6  ########################################
#######################################################################################################

chr <- "X"
bpstart <- 115200000
bpend <- 115700000

# agtr2_eQTL_files <- c("fastQTL_runs/Pancreas_FastQTL_X_114500000_116600000_no_age_fixed_slopeSE_eQTLs_mod.txt.gz"
#                         ,"Small_Intestine_Terminal_Ileum.v7.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
#                         ,"Stomach.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
#                         ,"Colon_Transverse.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
#                         ,"Colon_Sigmoid.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
# )
agtr2_eQTL_files <- c("Pancreas.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Small_Intestine_Terminal_Ileum.v7.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                      , "Lung.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Stomach.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Colon_Transverse.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
                        ,"Colon_Sigmoid.allpairs_fixed_chrX_114500000_116600000_eQTLs_mod.txt.gz"
)

tissues <- c("Pancreas", "Small_intestine", "Lung", "Stomach", "Colon_Transverse", "Colon_Sigmoid")

pdf("v7_plots/AGTR2_vs_eQTLs_overall_GTExV7.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1,2.1))
### AGTR2 overall ####
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs3788766_slc6a14.ld"
lead_snpname <- "chr23:115566839"
extra_y_height <- 5
ymax <- max(-log10(min(All_GWAS_pvalues$p))) + extra_y_height
ymin <- (-1)*(2/(get_eQTL_yrange2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000180772.6", chr, bpstart, bpend)+1))*ymax
plot_locuszoom(All_GWAS_pvalues, chr,bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
               r2colors=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"))
#plot_r2_legend("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "black")
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7879546_agtr2.ld"
lead_snpname <- as.character(lung_assoc_overall_UNC$MarkerName[which(lung_assoc_overall_UNC$p == min(lung_assoc_overall_UNC$p))])
lead_snpname <- gsub("chrX", "chr23", lead_snpname)
plot_locuszoom(lung_assoc_overall_UNC, chr,bpstart,bpend, yrange=c(-2,6),ld_filename, lead_snpname, transparency2, superimpose = T, yaxis="left",
               r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))
#plot_r2_legend("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","black")
plot_eQTL2(paste0(gtex_dir, agtr2_eQTL_files), "ENSG00000180772.6", chr, bpstart,bpend,window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir, agtr2_eQTL_files), "ENSG00000180772.6", chr, bpstart, bpend), 
           superimpose=T, yaxis="right",
           draw_dots = F, draw_linefit = T, title="AGTR2",lty=seq(1:length(tissues)), col=cols)
legend("top", tissues, title="AGTR2 Gene Expression eQTL P-values",
       lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')

# legend("topleft", title=expression("Lung GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"), bty='n')


dumbo = 0

agtr2_exon_starts <- c(115301957,115302220,115303498)/10^6
agtr2_exon_ends <- c(115302069,115302280,115306225)/10^6
draw_gene(y=-1.5+dumbo, agtr2_exon_starts,agtr2_exon_ends, "AGTR2 -->",bpstart/10^6,bpend/10^6)

slc6a14_exon_starts <- c(115567746,115568957,115572133,115573854,115574810,115576085,115577906,115582606,115584181,115585489,115586142,115586522,115588774,115589974)/10^6
slc6a14_exon_ends <- c(115567925,115569123,115572265,115574016,115574958,115576218,115578047,115582835,115584307,115585608,115586242,115586632,115588942,115592625)/10^6
draw_gene(y=-1.25+dumbo,slc6a14_exon_starts,slc6a14_exon_ends,"SLC6A14 -->",bpstart/10^6,bpend/10^6)

ct83_exon_starts <- c(115592852,115593942)/10^6
ct83_exon_ends <- c(115593174,115594194)/10^6
ct83_length <- length(ct83_exon_starts)
draw_gene(y=-1.85+dumbo,ct83_exon_starts,ct83_exon_ends,"<-- CT83",bpstart/10^6,bpend/10^6)


dev.off()



#######################################################################################################
#################################  SLC6A14 Males vs Females  ##########################################
#######################################################################################################


chr <- "X"
bpstart <- 115200000
bpend <- 115700000

slc6a14_eQTL_files <- c("fastQTL_runs/Pancreas_FastQTL_X_114500000_116600000_males_fixed_slopeSE_eQTLs_mod.txt.gz"
                        ,"fastQTL_runs/Pancreas_FastQTL_X_114500000_116600000_females_fixed_slopeSE_eQTLs_mod.txt.gz"
                        ,"HNE_slc6a14_eQTLs_males_mod.txt"
                        ,"HNE_slc6a14_eQTLs_females_mod.txt"
)
tissues <- c("Pancreas males (N=130)", "Pancreas females (N=90)", "HNE males (N=36)", "HNE females (N=40)")

pdf("v7_plots/SLC6A14_vs_eQTLs_male_female_eQTL_GTExV7_plus_76HNE.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+0,2.1+2)) # par()$mar
#par(mar=c(5.1,4.1,15.1,2.1)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs3788766_slc6a14.ld"
lead_snpname <- "chr23:115566839"
extra_y_height <- 5

ymax <- max(-log10(min(All_GWAS_pvalues$p))) + extra_y_height
#ymax=16
#ymin <- (-1)*(2/(get_eQTL_yrange("SLC6A14",tissues)+1))*ymax
ymin <- -4
# plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
#                r2colors=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"))
plot_locuszoom(All_GWAS_pvalues, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left",
               r2colors=c("purple","red","orange","green","cyan","blue","black"))

#plot_r2_legend("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "black")
# ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs7879546_agtr2.ld"
# lead_snpname <- as.character(lung_assoc_overall_UNC$MarkerName[which(lung_assoc_overall_UNC$p == min(lung_assoc_overall_UNC$p))])
# lead_snpname <- gsub("chrX", "chr23", lead_snpname)
# plot_locuszoom(lung_assoc_overall_UNC, "X",bpstart,bpend,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = T, yaxis="left",
#                r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))
#plot_r2_legend("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","black")
# plot_eQTL(genename="SLC6A14", tissues,chr="X",start=bpstart,end=bpend,window_size=25,ymax=get_eQTL_yrange("SLC6A14",tissues)+1, superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="SLC6A14",lty=seq(1:length(tissues)), colour=cols)
# legend("top", tissues, title="SLC6A14 eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')
plot_eQTL2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7",chr="X",start=bpstart,end=bpend,window_size=25,
           ymax=get_eQTL_yrange2(paste0(gtex_dir, slc6a14_eQTL_files), "ENSG00000087916.7", chr, bpstart, bpend)+1, 
           superimpose=T, yaxis="right",
           draw_dots = F, draw_linefit = T, title="SLC6A14",lty=seq(1:length(tissues)), colour=cols)
legend("top", tissues, title="SLC6A14 eQTL P-values",
       lty=seq(1:length(tissues)), lwd=rep(2,length(tissues)),col=cols, bty='n')


# legend("topleft", title=expression("Lung GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","#FF0000","#FF6600","#FFCC00","#CC9933","#CCCC33", "grey"), bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

# genename="slc6a14"
# gene_tissue_eqtls <- read.table(paste0(gtex_dir,i,"_",tolower(genename), "_eQTLs_mod.txt"), header=T, stringsAsFactors = F)
# ymax=get_eQTL_yrange("SLC6A14",tissues)+1
# pancreas_eqtl <- data.frame(snp=as.character(gene_tissue_eqtls$MarkerName), chr="X", pos=as.integer(gsub("X:(\\d*?)","\\1",gene_tissue_eqtls$MarkerName)),
#                             p=gene_tissue_eqtls$pval_nominal, MarkerName=gene_tissue_eqtls$MarkerName, stringsAsFactors = F)
# plot_locuszoom(pancreas_eqtl, "X", 115400000,115700000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency2, superimpose = F, yaxis="left",
#                r2colors=c("purple","#00FF33","#00CC00","#339933","#6666CC","#6699CC","grey"))


dumbo = 0

agtr2_exon_starts <- c(115301957,115302220,115303498)/10^6
agtr2_exon_ends <- c(115302069,115302280,115306225)/10^6
draw_gene(y=-1.5+dumbo, agtr2_exon_starts,agtr2_exon_ends, "AGTR2 -->",bpstart/10^6,bpend/10^6)

slc6a14_exon_starts <- c(115567746,115568957,115572133,115573854,115574810,115576085,115577906,115582606,115584181,115585489,115586142,115586522,115588774,115589974)/10^6
slc6a14_exon_ends <- c(115567925,115569123,115572265,115574016,115574958,115576218,115578047,115582835,115584307,115585608,115586242,115586632,115588942,115592625)/10^6
draw_gene(y=-1.25+dumbo,slc6a14_exon_starts,slc6a14_exon_ends,"SLC6A14 -->",bpstart/10^6,bpend/10^6)

ct83_exon_starts <- c(115592852,115593942)/10^6
ct83_exon_ends <- c(115593174,115594194)/10^6
ct83_length <- length(ct83_exon_starts)
draw_gene(y=-1.85+dumbo,ct83_exon_starts,ct83_exon_ends,"<-- CT83",bpstart/10^6,bpend/10^6)

#dev.off()



####### Where are the enhancers and promoters? ########
### None to be found in any of the tissues for chrX


ystart = 8
step = 0.3
trackheight = 0.4
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=23, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight,
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=23, coord_range=c(bpstart,bpend), ystart+trackheight, trackheight,
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(bpend/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()




#######################################################################################################
############################################  PRSS1  ##################################################
#######################################################################################################

chr <- "7"
bpstart <- 142350000
bpend <- 142600000

pdf("PRSS1_vs_eQTLs_Rplot8.pdf")
#pdf("PRSS1_locuszoom.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+2,2.1+6)) # par()$mar
#par(mar=c(5.1,4.1,4.1,2.1)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs10273639_prss1.ld"
lead_snpname <- "chr7:142456928"
extra_y_height <- 2
region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 7 & All_GWAS_pvalues$pos >= 142350000 & All_GWAS_pvalues$pos <= 142600000)
#ymax <- max(-log10(min(region$p))) + extra_y_height
ymax <- 15
#ymin <- (-1)*(2/(get_eQTL_yrange("PRSS1",tissues)+1))*ymax
ymin <- (-1)*(2/(3+1))*ymax
plot_locuszoom(All_GWAS_pvalues, "7",142350000,142600000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")
# plot_eQTL("PRSS1", tissues,"7",142350000,142600000,window_size=25,ymax=get_eQTL_yrange("PRSS1",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="PRSS1",lty=seq(1:length(tissues)), colour=cols)
plot_eQTL("PRSS1", tissues,"7",142350000,142600000,window_size=25,ymax=3, superimpose=T, yaxis="right",
          draw_dots = F, draw_linefit = T, title="PRSS1",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="PRSS1 Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

MTRNR2L6_exon_starts <- c(142374130)/10^6
MTRNR2L6_exon_ends <- c(142375525)/10^6
PRSS1_exon_starts <- c(142457318,142458405,142459624,142460281,142460718)/10^6
PRSS1_exon_ends <- c(142457375,142458565,142459878,142460418,142460927)/10^6
PRSS3P2_exon_starts <- c(142478756,142479908,142481126,142481775,142482211)/10^6
PRSS3P2_exon_ends <- c(142478879,142480068,142481380,142481912,142482399)/10^6
EPHB6_exon_starts <- c(142552775,142558809,142559774,142560513,142560884,142561723,142563229,142563714,142564235,142564660,142565362,142565751,142565995,142566246,142566726,142567569,142567966,142568282,142568548)/10^6
EPHB6_exon_ends <- c(142553147,142558944,142559891,142560591,142561085,142562504,142563385,142564071,142564360,142564823,142565477,142565804,142566115,142566494,142566900,142567719,142568160,142568438,142568847)/10^6
TRPV6_exon_starts <- c(142568955,142570124,142571200,142571828,142572243,142572656,142572830,142573220,142573510,142574160,142574491,142574894,142575403,142575681,142583133)/10^6
TRPV6_exon_ends <- c(142569742,142570231,142571469,142571895,142572409,142572733,142572917,142573433,142573657,142574336,142574590,142575032,142575526,142575779,142583490)/10^6

draw_gene(y=-0.75, MTRNR2L6_exon_starts, MTRNR2L6_exon_ends, "MTRNR2L6 -->",bpstart/10^6, bpend/10^6)
draw_gene(y=-0.75, PRSS1_exon_starts, PRSS1_exon_ends, "PRSS1 -->",bpstart/10^6, bpend/10^6)
draw_gene(y=-1.5, PRSS3P2_exon_starts, PRSS3P2_exon_ends, "PRSS3P2 -->",bpstart/10^6, bpend/10^6)
draw_gene(y=-0.75, EPHB6_exon_starts, EPHB6_exon_ends, "EPHB6 -->",bpstart/10^6, bpend/10^6)
draw_gene(y=-1.5, TRPV6_exon_starts, TRPV6_exon_ends, "<-- TRPV6",bpstart/10^6, bpend/10^6)

# dev.off()




####### Where are the enhancers and promoters? ########

ystart = 3
step = 0.35
trackheight = 0.3
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=7, coord_range=c(142350000,142600000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=7, coord_range=c(142350000,142600000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(142610000/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()



#######################################################################################################
##############################################  CEBPB  ################################################
#######################################################################################################


pdf("CEBPB_vs_eQTLs_Rplot4.pdf")
#pdf("UBE2V1_vs_eQTLs_Rplot4.pdf")
#pdf("TMEM189_vs_eQTLs_Rplot4.pdf")
#pdf("RP11-112L6.2_vs_eQTLs_Rplot3.pdf")
#pdf("RP11-112L6.3_vs_eQTLs_Rplot4.pdf")
#pdf("RP11-112L6.4_vs_eQTLs_Rplot4.pdf")
#pdf("RP11-290F20.3_vs_eQTLs_Rplot4.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+14,2.1+8)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs867489_cebpb.ld"
lead_snpname <- "chr20:48833957"
extra_y_height <- 5
region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 20 & All_GWAS_pvalues$pos >= 48700000 & All_GWAS_pvalues$pos <= 48900000)
ymax <- max(-log10(min(region$p))) + extra_y_height
ymin <- (-1)*(2/(get_eQTL_yrange("CEBPB",tissues)+1))*ymax
plot_locuszoom(All_GWAS_pvalues, "20",48700000,48900000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")
plot_eQTL("CEBPB", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("CEBPB",tissues), superimpose=T, yaxis="right",
          draw_dots = F, draw_linefit = T, title="CEBPB",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="CEBPB Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

UBE2V1_exon_starts <- c(48697660,48700665,48713019,48713208,48732021)/10^6
UBE2V1_exon_ends <- c(48699451,48700791,48713071,48713357,48732496)/10^6
TMEM189_UBE2V1_exon_starts <- c(48697660,48700665,48713208,48744511,48746082,48747392,48760038,48770053)/10^6
TMEM189_UBE2V1_exon_ends <- c(48699451,48700791,48713357,48744724,48746227,48747484,48760158,48770335)/10^6
TMEM189_exon_starts <- c(48740273,48744511,48746082,48747392,48760038,48770053)/10^6
TMEM189_exon_ends <- c(48741716,48744724,48746227,48747484,48760158,48770335)/10^6
LINC01273_exon_starts <- c(48789106,48789459,48792817)/10^6
LINC01273_exon_ends <- c(48789325,48789554,48793208)/10^6
CEBPB_AS1_exon_starts <- c(48801139,48803549,48805079,48805992,48806534,48807804,48808200)/10^6
CEBPB_AS1_exon_ends <- c(48802957,48803664,48805164,48806098,48806822,48807959,48808606)/10^6
CEBPB_exon_starts <- c(48807119)/10^6
CEBPB_exon_ends <- c(48809227)/10^6
LINC01272_exon_starts <- c(48884022,48892182,48894027,48894713)/10^6
LINC01272_exon_ends <- c(48884200,48892272,48894150,48896332)/10^6

draw_gene(y=-2/3, UBE2V1_exon_starts, UBE2V1_exon_ends, "<-- UBE2V1", cex=0.5)
draw_gene(y=(-2/3)*2, TMEM189_UBE2V1_exon_starts, TMEM189_UBE2V1_exon_ends, "", cex=0.5)
text(UBE2V1_exon_ends[length(UBE2V1_exon_ends)], -0.75, "<-- TMEM189-UBE2V1", pos=4, cex=0.5, col="red")
draw_gene(y=(-2/3)*3, TMEM189_exon_starts, TMEM189_exon_ends, "", cex=0.5)
text(TMEM189_exon_ends[length(TMEM189_exon_ends)],-2.1,"<-- TMEM189",pos=2,cex=0.5,col="red")
draw_gene(y=-2/3, CEBPB_exon_starts, CEBPB_exon_ends, "CEBPB -->", cex=0.5)
draw_gene(y=(-2/3)*2, CEBPB_AS1_exon_starts, CEBPB_AS1_exon_ends, "<-- CEBPB-AS1", cex=0.5)
draw_gene(y=-2, LINC01273_exon_starts, LINC01273_exon_ends, "LINC01273 -->", cex=0.5)
draw_gene(y=-2/3, LINC01272_exon_starts, LINC01272_exon_ends, "", cex=0.5)
text(LINC01272_exon_ends[length(LINC01272_exon_ends)], -0.8, "LINC01272 -->", pos=2, cex=0.5, col="red")

#dev.off()



######### eQTLs for UBE2V1 drawn instead ###########
# plot_eQTL("UBE2V1", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("UBE2V1",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="UBE2V1",lty=seq(1:length(tissues)), colour=cols)


######### eQTLs for TMEM189 drawn instead ###########
# plot_eQTL("TMEM189", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("TMEM189",tissues), superimpose=T, yaxis="right",
#             draw_dots = F, draw_linefit = T, title="TMEM189",lty=seq(1:length(tissues)), colour=cols)

######### eQTLs for RP11-112L6.2 drawn instead ###########
# No expression in any of the tissues #
# plot_eQTL("RP11-112L6.2", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("RP11-112L6.2",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="RP11-112L6.2",lty=seq(1:length(tissues)), colour=cols)


######### eQTLs for RP11-112L6.3 drawn instead ###########
# plot_eQTL("RP11-112L6.3", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("RP11-112L6.3",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="RP11-112L6.3",lty=seq(1:length(tissues)), colour=cols)


######### eQTLs for RP11-112L6.4 drawn instead ###########
# plot_eQTL("RP11-112L6.4", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("RP11-112L6.4",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="RP11-112L6.4",lty=seq(1:length(tissues)), colour=cols)

# ######### eQTLs for RP11-290F20.3 drawn instead ###########
# plot_eQTL("RP11-290F20.3", tissues,"20",48700000,48900000,window_size=25,ymax=get_eQTL_yrange("RP11-290F20.3",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="RP11-112L6.4",lty=seq(1:length(tissues)), colour=cols)



####### Where are the enhancers and promoters? ########

ystart = get_eQTL_yrange("CEBPB",tissues) + 1
step = 0.5
trackheight = 0.45
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=20, coord_range=c(48700000,48900000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=20, coord_range=c(48700000,48900000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(48910000/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()




#######################################################################################################
#############################################  CFTR  ##################################################
#######################################################################################################

#pdf("CFTR_vs_eQTLs_Rplot3.pdf")
pdf("CFTR_vs_eQTLs_Rplot4.1.pdf")
par(oma=c(0,0,0,2))
par(mar=c(5.1,4.1,4.1+14,2.1+7)) # par()$mar
ld_filename <- "/Users/naim/Downloads/locuszoom/data/1000G/genotypes/2012-03/EUR/rs213958_MI_CFTR_SNP_r2.ld"
lead_snpname <- "chr7:117223764"
extra_y_height <- 2
# region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 7 & All_GWAS_pvalues$pos >= 116106355 & All_GWAS_pvalues$pos <= 118105720)
region <- subset(All_GWAS_pvalues, All_GWAS_pvalues$chr %in% 7 & All_GWAS_pvalues$pos >= 116850000 & All_GWAS_pvalues$pos <= 117500000)
ymax <- max(-log10(min(region$p))) + extra_y_height
ymin <- (-1)*(2/(get_eQTL_yrange("CFTR",tissues)+1))*ymax
# plot_locuszoom(All_GWAS_pvalues, "7",116106355,118105720,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")
plot_locuszoom(All_GWAS_pvalues, "7",116850000,117500000,yrange=c(ymin,ymax),ld_filename, lead_snpname, transparency, superimpose = F, yaxis="left")
# plot_eQTL("CFTR", tissues,"7",116106355,118105720,window_size=25,ymax=get_eQTL_yrange("CFTR",tissues), superimpose=T, yaxis="right",
#           draw_dots = F, draw_linefit = T, title="CFTR",lty=seq(1:length(tissues)), colour=cols)
plot_eQTL("CFTR", tissues,"7",116850000,117500000,window_size=25,ymax=get_eQTL_yrange("CFTR",tissues), superimpose=T, yaxis="right",
          draw_dots = F, draw_linefit = T, title="CFTR",lty=seq(1:length(tissues)), colour=cols)
# legend("topleft", tissues, title="CFTR Gene Expression eQTL P-values",
#        lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')
# legend("topright", title=expression("MI GWAS" ~ r^{2}), legend=c("REF SNP", "0.8", "0.6", "0.4", "0.2","<0.2", "Unknown"),
#        pch=rep(16,7), col=c("purple","red","orange","green","cyan","blue","black"), bty='n')

CAV2_exon_starts <- c(116139654,116140313,116146024)/10^6
CAV2_exon_ends <- c(116139985,116140501,116148595)/10^6
CAV1_exon_starts <- c(116164838,116166578,116198999)/10^6
CAV1_exon_ends <- c(116165146,116166743,116201239)/10^6
LINC01510_exon_starts <- c(116203647,116211171,116213494,116240844,116254652)/10^6
LINC01510_exon_ends <- c(116204222,116211284,116213798,116240988,116254874)/10^6
MET_exon_starts <- c(116312412,116371721,116380003,116380905,116395408,116397490,116397691,116398512,116399444,116403103,116409698,116411551,116411902,116414934,116417442,116418829,116422041,116423357,116435708,116435940)/10^6
MET_exon_ends <- c(116312631,116371913,116380138,116381079,116395569,116397593,116397828,116398674,116399544,116403322,116409845,116411708,116412043,116415165,116417523,116419011,116422151,116423523,116435845,116438440)/10^6
CAPZA2_exon_starts <- c(116502562,116528180,116533047,116538825,116544230,116546316,116550286,116552122,116556113,116557780)/10^6
CAPZA2_exon_ends <- c(116502704,116528244,116533099,116538889,116544437,116546396,116550365,116552194,116556176,116559313)/10^6
ST7_AS1_exon_starts <- c(116592500)/10^6
ST7_AS1_exon_ends <- c(116594388)/10^6
ST7_exon_starts <- c(116593380,116739815,116759614,116769846,116770544,116771938,116776134,116778488,116810915,116829374,116830887,116849840,116859137,116861976,116869815)/10^6
ST7_exon_ends <- c(116593745,116739898,116759774,116769901,116770660,116772014,116776289,116778586,116811030,116829447,116830990,116849991,116859230,116862116,116870075)/10^6
ST7_OT4_exon_starts <- c(116593952,116594514,116595027,116596476,116599206)/10^6
ST7_OT4_exon_ends <- c(116594163,116594733,116595207,116596777,116599867)/10^6
MIR6132_exon_starts <- c(116660264)/10^6
MIR6132_exon_ends <- c(116660373)/10^6
ST7_AS2_exon_starts <- c(116712125,116712507,116713873,116720918,116730775)/10^6
ST7_AS2_exon_ends <- c(116712415,116713194,116713980,116721126,116730835)/10^6
ST7_AS2_exon_starts2 <- c(116752345,116757954,116768350,116770581,116784220,116785449)/10^6
ST7_AS2_exon_ends2 <- c(116752721,116758037,116768397,116770741,116784307,116785646)/10^6
ST7_OT3_exon_starts <- c(116822734,116824178,116829374,116830186,116830887,116831688,116838373,116839125,116849840)/10^6
ST7_OT3_exon_ends <- c(116822799,116824253,116829447,116830322,116830990,116831815,116838594,116839198,116849991)/10^6
WNT2_exon_starts <- c(116916685,116937665,116955124,116960620,116962960)/10^6
WNT2_exon_ends <- c(116918438,116937930,116955402,116960847,116963343)/10^6
ASZ1_exon_starts <- c(117003275,117007405,117008665,117019991,117021064,117022122,117023039,117024779,117025751,117060216,117062290,117066889,117067409)/10^6
ASZ1_exon_ends <- c(117003802,117007492,117008771,117020101,117021121,117022198,117023164,117024914,117025863,117060328,117062413,117066989,117067577)/10^6
CFTR_exon_starts <- c(117120016,117144306,117149087,117170952,117174329,117175301,117176601,117180153,117182069,117188694,117199517,117227792,117230406,117231987,117234983,117242879,117243585,117246727,117250572,117251634,117254666,117267575,117282491,117292895,117304741,117305512,117306961)/10^6
CFTR_exon_ends <- c(117120201,117144417,117149196,117171168,117174419,117175465,117176727,117180400,117182162,117188877,117199709,117227887,117230493,117232711,117235112,117242917,117243836,117246807,117250723,117251862,117254767,117267824,117282647,117292985,117304914,117305618,117308718)/10^6
CTTNB2_exon_starts <- c(117350705,117358071,117359557,117361120,117364600,117365105,117368142,117374966,117375322,117385884,117386066,117396608,117397928,117400488,117407112,117417564,117420494,117422915,117424304,117431181,117450818,117501262,117513388)/10^6
CTTNB2_exon_ends <- c(117351836,117358173,117359690,117361184,117364786,117365311,117368321,117375154,117375475,117385984,117386153,117396688,117398024,117400764,117407230,117417819,117420645,117423015,117424508,117432835,117451043,117501370,117513561)/10^6
LSM8_exon_starts <- c(117824085,117825707,117828331,117831965)/10^6
LSM8_exon_ends <- c(117824308,117825748,117828459,117844093)/10^6
ANKRD7_exon_starts <- c(117864711,117874484,117874754,117876094,117876843,117879962,117882402)/10^6
ANKRD7_exon_ends <- c(117865063,117874599,117874928,117876201,117876980,117880052,117882784)/10^6

# draw_gene(y=-2, CAV2_exon_starts, CAV2_exon_ends, "CAV2",cex=0.5)
# draw_gene(y=-1.5, CAV1_exon_starts, CAV1_exon_ends, "CAV1 -->",cex=0.5)
# draw_gene(y=-1, LINC01510_exon_starts, LINC01510_exon_ends, "LINC01510", cex=0.5)
# draw_gene(y=-2, MET_exon_starts,MET_exon_ends,"MET",cex=0.5)
# draw_gene(y=-1.5, CAPZA2_exon_starts, CAPZA2_exon_ends, "CAPZA2",cex=0.5)
# draw_gene(y=-1, ST7_AS1_exon_starts, ST7_AS1_exon_ends,"ST7-AS1",cex=0.5)
# draw_gene(y=-0.5, ST7_exon_starts, ST7_exon_ends, "ST7", cex=0.5)
draw_gene(y=-1, WNT2_exon_starts,WNT2_exon_ends,"WNT2",cex=0.5)
# draw_gene(y=-2, ST7_OT4_exon_starts, ST7_OT4_exon_ends, "ST7-OT4",cex=0.5)
# draw_gene(y=-2.2, MIR6132_exon_starts, MIR6132_exon_ends, "MIR6132",cex=0.5)
# draw_gene(y=-1.8, ST7_AS2_exon_starts,ST7_AS2_exon_ends, "ST7-AS2",cex=0.5)
# draw_gene(y=-0.75, ST7_AS2_exon_starts2, ST7_AS2_exon_ends2, "ST7-AS2",cex=0.5)
draw_gene(y=-1.5,ASZ1_exon_starts,ASZ1_exon_ends,"ASZ1",cex=0.5)
draw_gene(y=-0.5, CFTR_exon_starts,CFTR_exon_ends,"CFTR",cex=0.5)
draw_gene(y=-1,CTTNB2_exon_starts,CTTNB2_exon_ends, "CTTNB2",cex=0.5)
# draw_gene(y=-1, LSM8_exon_starts,LSM8_exon_ends,"LSM8",cex=0.5)
# draw_gene(y=-1.5,ANKRD7_exon_starts,ANKRD7_exon_ends,"ANKRD7",cex=0.5)

# dev.off()




####### Where are the enhancers and promoters? ########

ystart = 5.2
step = 0.5
trackheight = 0.4
#plot(x=c(205.86,205.94),y=c(ystart,ystart-step*nrow(relevant_tissues)),type="n")
counter=0
for(i in relevant_tissues[,1]) {
  counter=counter+1
  enh_index <- which(colnames(enhancers) %in% i)
  enh_tissue_i <- which(enhancers[,enh_index] != 0)
  colour_regions_track(chr=7, coord_range=c(116850000,117500000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(enhancers)[enh_tissue_i], positive_regions_col="darkgoldenrod", units="Mb")
  prom_index <- which(colnames(promoters) %in% i)
  prom_tissue_i <- which(promoters[,prom_index] != 0)
  colour_regions_track(chr=7, coord_range=c(116850000,117500000), ystart+trackheight, trackheight, 
                       positive_regions=rownames(promoters)[prom_tissue_i], positive_regions_col="red", units="Mb")
  text(117510000/10^6, ystart + trackheight/2, labels=relevant_tissues[counter,2], pos=4)
  ystart <- ystart + step
}

dev.off()






########## Separate legend ############
pdf("eQTL_legend.pdf")
par(oma=c(0,0,0,0))
par(mar=c(5.1,4.1,4.1,2.1)) # par()$mar
plot(c(1,1),c(1,1),type='n', xaxt='n',yaxt='n', xlab='', ylab='')
legend("center", tissues, title="Gene Expression eQTL P-values",
       lty=seq(1:length(tissues)), lwd=c(2,2,2,2,2),col=cols, bty='n')
dev.off()
