#!/usr/bin/R
# This version simply runs PC-AiR and saves the resulting eigenvalues and eigenvectors

args=(commandArgs(TRUE))


bplink_filename <- as.character(args[1]) # Binary plink filename w/out the .bed/.bim/.fam extension
kinship_filename <- as.character(args[2])
out_prefix <- as.character(args[3])
# reference_panels_included <- as.logical(args[4])  # feature to be implemented in later versions

# source("https://bioconductor.org/biocLite.R")
# biocLite("GWASTools")
# biocLite("GENESIS")
# biocLite("SNPRelate")


library(GWASTools)
library(GENESIS)
library(SNPRelate)

#unlink("genotype.gds",force=T)


print("Converting binary PLINK file to GDS format")
snpgdsBED2GDS(bed.fn = paste0(bplink_filename,".bed"),
              bim.fn = paste0(bplink_filename,".bim"),
              fam.fn = paste0(bplink_filename,".fam"),
              out.gdsfn = paste0(bplink_filename,".gds"))
print("Reading genotype data")
geno <- GdsGenotypeReader(filename = paste0(bplink_filename,".gds"))
genoData <- GenotypeData(geno)

iids <- getScanID(genoData)
print("Reading KING kinship coefficients")
KINGmat <- king2mat(file.kin0 = paste0(kinship_filename,".kin0"),
                    file.kin = paste0(kinship_filename,".kin"), iids = iids)
# or run KING using snpgdsIBDKING() function

print("Running PC-AiR")
mypcair <- pcair(genoData = genoData, kinMat = KINGmat, divMat = KINGmat)

print("Generating PC-AiR plots for the first 3 PCs")
pdf(paste0(out_prefix,"_general_pca_plot_1of3.pdf"))
plot(mypcair, vx = 1, vy = 2)
dev.off()
pdf(paste0(out_prefix,"_general_pca_plot_2of3.pdf"))
plot(mypcair, vx = 2, vy = 3)
dev.off()
pdf(paste0(out_prefix,"_general_pca_plot_3of3.pdf"))
plot(mypcair, vx = 3, vy = 4)
dev.off()


eigenvalues <- mypcair$values
eigenvectors <- mypcair$vectors


str(mypcair)


print("Saving the eigenvalues and eigenvectors")
# Save eigenvalues and eigenvectors:
write.table(eigenvalues, paste0(out_prefix,"_eigenvalues.txt"), quote=F, row.names=T, col.names=F)
write.table(eigenvectors, paste0(out_prefix,"_eigenvectors.txt"), quote=F, row.names=T, col.names=F)
write.table(mypcair$rels, paste0(out_prefix,"_relateds.txt"), quote=F, row.names=F,col.names=F)
write.table(mypcair$unrels, paste0(out_prefix,"_unrelateds.txt"),quote=F,row.names=F,col.names=F)


print("Done")

# In my later R versions, will try to enable the plotting of cases and controls, and if included, the reference panels (v3)





#====================================================================
# FUNCTIONS
#====================================================================
boxArea <- function(dataPoints) {
  # PRE: given 3 vectors of data points like 3 PCs in data.frame format
  # POST: returns the xrange, yrange and zrange of the 3D box as a vector 
  #     in the following format: c(xmin, xmax, ymin, ymax, zmin, zmax)
  xmin <- min(dataPoints[,1])
  xmax <- max(dataPoints[,1])
  ymin <- min(dataPoints[,2])
  ymax <- max(dataPoints[,2])
  zmin <- min(dataPoints[,3])
  zmax <- max(dataPoints[,3])
  return(c(xmin, xmax, ymin, ymax, zmin, zmax))
}

isInRange <- function(range, point) (point >= range[1] & point <= range[2])


isInBox <- function(boxPoints, samplePoints) {
  # PRE: give the 3 PCs from hapmap data for specific ethnic group in boxPoints and 
  #       the sample individuals to be tested for belonging in the
  #       cluster under samplePoints
  # POST: returns boolean vector on whether individual is within ethnic group
  inBox <- logical(length=dim(samplePoints)[1])
  box <- boxArea(boxPoints)
  xrange <- c(box[1], box[2])
  yrange <- c(box[3], box[4])
  zrange <- c(box[5], box[6])
  for (i in 1:dim(samplePoints)[1]) {
    if (isInRange(xrange, samplePoints[i,1]) & 
        isInRange(yrange, samplePoints[i,2]) &
        isInRange(zrange, samplePoints[i,3]))
    {
      inBox[i] <- TRUE
    } else {
      inBox[i] <- FALSE
    }
  }
  return(inBox)
}

getOutliers <- function(eigen_set, eigen_sd_set, fid_column, iid_column, eigenvalue_columns) {
  final_outliers_iid <- NULL
  for(i in eigenvalue_columns) {
    left_limit <- mean(as.numeric(as.character(eigen_sd_set[,i]))) - 6*as.numeric(as.character(sd(eigen_sd_set[,i])))
    right_limit <- mean(as.numeric(as.character(eigen_sd_set[,i]))) + 6*as.numeric(as.character(sd(eigen_sd_set[,i])))
    outliers_vec_i <- subset(eigen_set, as.numeric(as.character(eigen_set[,i])) < left_limit | as.numeric(as.character(eigen_set[,i])) > right_limit)
    
    final_outliers_iid <- unique(c(final_outliers_iid, as.character(outliers_vec_i[,iid_column])))
  }
  final_outliers_fid <- getFID(final_outliers_iid)
  outliers <- data.frame(FID=final_outliers_fid, IID=final_outliers_iid)
  return(outliers)
}

isInBox_Nsd <- function(boxPoints, samplePoints, N) {
  # PRE: give the 3 PCs from hapmap data for specific ethnic group in boxPoints and 
  #       the sample individuals to be tested for belonging in the
  #       cluster under samplePoints
  # POST: returns boolean vector on whether individual is within ethnic group +/- N sd
  inBox <- logical(length=dim(samplePoints)[1])
  box <- boxArea(boxPoints)
  xrange <- c(box[1] - N*sd(boxPoints[,1]), box[2] + N*sd(boxPoints[,1]))
  yrange <- c(box[3] - N*sd(boxPoints[,2]), box[4] + N*sd(boxPoints[,2]))
  zrange <- c(box[5] - N*sd(boxPoints[,3]), box[6] + N*sd(boxPoints[,3]))
  for (i in 1:dim(samplePoints)[1]) {
    if (isInRange(xrange, samplePoints[i,1]) & 
        isInRange(yrange, samplePoints[i,2]) &
        isInRange(zrange, samplePoints[i,3]))
    {
      inBox[i] <- TRUE
    } else {
      inBox[i] <- FALSE
    }
  }
  return(inBox)
}

getFID <- function(IIDs) {
  fid <- character(length(IIDs))
  for (i in 1:length(IIDs)) {
    fid[i] <- as.character(subset(peoples$FID, peoples$IID %in% as.character(IIDs[i])))
  }
  return(fid)
}
#====================================================================
# END OF FUNCTIONS
#====================================================================


if(reference_panels_included) {
  the1kgpop <- read.table("~/1000Genomes_phase3_for_BEAGLE_v4/integrated_call_samples_v3.20130502.ALL.panel", header=T, stringsAsFactors = F)
  the1kgped <- read.table("~/1000Genomes_phase3_for_BEAGLE_v4/integrated_call_samples.20130502.ALL.ped", sep="\t", header=T, stringsAsFactors = F)
  populations <- levels(as.factor(the1kgpop$super_pop))
  colours <- c("green", "purple", "black", "brown", "blue")
  samples <- read.table("23_UKomni25for702_dups_nameUpdate.fam", header=F, stringsAsFactors=F)
  the1kg <- read.table("chr.ALL.1kg.phase3.v5a.bplink.dups_rm.fam", stringsAsFactors=F)
  eigen_1kg <- subset(eigenvectors, rownames(eigenvectors) %in% the1kg[,2])
  eigen_sample <- subset(eigenvectors, rownames(eigenvectors) %in% samples[,2])
  samples_pop <- cbind(sample=samples[,2], pop="Sample", super_pop="Sample", gender=ifelse(samples[,5]==1,"male","female"))
  allpops <- rbind(the1kgpop, samples_pop)
  populations <- c(populations, "Sample")
  colours <- c(colours, "red")
 
  # Hapmap version:
  # hapmappop <- read.table("~/hapmap_data/hapmap3r3_b36_relationship.txt",header=T,stringsAsFactors=F)
  # populations <- levels(as.factor(hapmappop$population))
  # hapmap <- read.table("hapmap_hg19_related_rm_ceu_tsi.fam", stringsAsFactors = F)
  # eigen_hapmap <- subset(eigenvectors, rownames(eigenvectors) %in% hapmap[,2])
  # hapmappop2 <- cbind(sample=hapmappop$IID, pop=hapmappop$population, super_pop=hapmappop$population,
  #                    gender=ifelse(hapmappop$sex==1,"male","female"))
  # allpops <- rbind(hapmappop2, samples_pop)
  # populations <- c(populations, "Sample")
  # colours <- c("green", "purple", "black", "brown", "blue",
  #             "cadetblue", "pink", 54, 53, 55, "orange", "red")
  # 
  # populations <- c("CEU", "TSI", "Sample")
  # colours <- colours[c(2,10,12)]
    
  draw2DPCA <- function() {   
    par(mfrow=c(2,2))
    # Eigenvectors 1 vs 2
    axis1 <- 1
    axis2 <- 2
    plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
    for(i in 1:length(populations)) {
      pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
      eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
      points(eigen_i[,c(axis1,axis2)], col=colours[i])
    }
    #  legend("topleft", c(populations), col=colours, pch=21)
    
    
    # Eigenvectors 1 vs 3
    axis1 <- 1
    axis2 <- 3
    plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
    for(i in 1:length(populations)) {
      pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
      eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
      points(eigen_i[,c(axis1,axis2)], col=colours[i])
    }
    #  legend("topleft", c(populations), col=colours, pch=21)
    
    
    # Eigenvectors 2 vs 3
    axis1 <- 2
    axis2 <- 3
    plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
    for(i in 1:length(populations)) {
      pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
      eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
      points(eigen_i[,c(axis1,axis2)], col=colours[i])
    }
    #  legend("topleft", c(populations), col=colours, pch=21)
    
    # Draw legend
    plot(c(0,1),c(0,length(populations)),type='n',xlab='',ylab='')
    x1<-rep(0,length(populations))
    x2<-rep(0.2,length(populations))
    y1=0:(length(populations)-1)
    y2=1:(length(populations))
    rect(x1,y1,x2,y2,col=colours)
    for(i in 1:length(populations)) {
      text(0.4,(i-1)+0.5,populations[i])
    }
  }
  
  
  # Want to determine ethnicity of each individual
  pop_i_samples <- list()
  esample_individuals <- eigen_sample
  print(paste("Number of individuals per population"))
  for (i in 1:(length(populations)-1)){ # all 1KG populations except Samples
    pop_i <- subset(the1kgpop, the1kgpop$super_pop %in% populations[i])
    eigen_i <- subset(eigen_1kg, rownames(eigen_1kg) %in% pop_i$sample)
    individual_indices <- isInBox(eigen_i[,1:3], esample_individuals[,1:3])
    print(paste(populations[i],sum(individual_indices)))
    pop_i_samples[[i]] <- as.character(rownames(esample_individuals)[individual_indices])
    names(pop_i_samples)[i] <- populations[i]
    if (sum(individual_indices)>0) esample_individuals <- esample_individuals[-which(individual_indices),] # remove sample individuals to prevent duplicates
  }
  sapply(pop_i_samples,length)
  sum(sapply(pop_i_samples,length))
  
  # who has not been classified?
  rownames(esample_individuals)
  
  # include last category of unclassified individuals:
  pop_i_samples[[length(pop_i_samples)+1]] <- as.character(rownames(esample_individuals))
  names(pop_i_samples)[length(pop_i_samples)] <- "Unknown"

  # Assign ethnicities per sample:
  ethnicity <- NULL
  fid <- NULL
  iid <- NULL
  for(i in 1:length(pop_i_samples)) {
    ethnicity <- c(ethnicity,rep(names(pop_i_samples)[i], length(pop_i_samples[[i]])))
    iid <- c(iid, pop_i_samples[[i]])
    for(j in 1:length(pop_i_samples[[i]])) {
      fid <- c(fid, subset(samples[,1], samples[,2] %in% pop_i_samples[[i]][j])[1])
    }
  }
  ind_ethnicities <- data.frame(FID=fid, IID=iid, Ethnicity=ethnicity)
  write.table(ind_ethnicities, "ethnicities.txt",
              quote=F, row.names=F, col.names=T)
  
  
  # Redraw and mark where the Unknowns are:
  populations <- c(populations, "Unknown")
  unknowns <- subset(ind_ethnicities, ind_ethnicities$Ethnicity %in% "Unknown")
  allpops[which(allpops[,1] %in% unknowns$IID),2] <- "Unknown"
  allpops[which(allpops[,1] %in% unknowns$IID),3] <- "Unknown"
  colours <- c(colours, "orange")

}


# Individuals' list to do 2nd PCA round:
samples_for_2ndPCA <- subset(ind_ethnicities, ind_ethnicities$Ethnicity %in% c("AMR", "EUR", "Unknown"))[,1:2]
unrelated_1kg <- subset(the1kg, the1kg[,2] %in% rownames(eigen_1kg)) # unrelated 1KG set
the1kgpopEUR <- subset(the1kgpop, the1kgpop$super_pop %in% "EUR")
unrelated_1kg_EUR <- subset(unrelated_1kg, unrelated_1kg[,2] %in% the1kgpopEUR$sample)
unrelated_1kg_EUR <- unrelated_1kg_EUR[,1:2]
colnames(unrelated_1kg_EUR) <- colnames(samples_for_2ndPCA)
samples_for_2ndPCA <- rbind(samples_for_2ndPCA, unrelated_1kg_EUR)

write.table(samples_for_2ndPCA, "samples_for_2ndPCA.txt", quote=F, row.names=F, col.names=F)



####### version 2.0 stop here ####



# Identify outliers:
peoples <- rbind(hapmap[,1:2], samples[,1:2])
colnames(peoples) <- c("FID", "IID")
eigen_iids <- matrix(ncol=1, rownames(eigenvectors))
colnames(eigen_iids) <- "IID"
newFIDs <- NULL
for(i in 1:nrow(eigen_iids)) {
  newFIDs <- c(newFIDs, peoples[which(peoples[,'IID'] %in% eigen_iids[i,'IID']),'FID'])
}
ordered_eigenIDs <- cbind(newFIDs, eigen_iids)

eigenvectors2 <- cbind(ordered_eigenIDs, eigenvectors)

(ctrl_sd_outliers <- getOutliers(eigen_set=eigenvectors2, eigen_sd_set=eigen_hapmap, fid_column=1, iid_column=2, eigenvalue_columns=3:5))
eigen_outliers <- subset(eigen_sample, rownames(eigen_sample) %in% as.character(ctrl_sd_outliers[,'IID']))


######## MARK OUTLIERS IN PLOT #######


par(mfrow=c(2,2))
# Eigenvectors 1 vs 2
axis1 <- 1
axis2 <- 2
plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
for(i in 1:length(populations)) {
  pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
  eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
  points(eigen_i[,c(axis1,axis2)], col=colours[i])
}
points(eigen_outliers[,c(axis1,axis2)], col="turquoise", pch=19, cex=0.3)
#  legend("topleft", c(populations), col=colours, pch=21)


# Eigenvectors 1 vs 3
axis1 <- 1
axis2 <- 3
plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
for(i in 1:length(populations)) {
  pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
  eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
  points(eigen_i[,c(axis1,axis2)], col=colours[i])
}
points(eigen_outliers[,c(axis1,axis2)], col="turquoise", pch=19, cex=0.3)
#  legend("topleft", c(populations), col=colours, pch=21)


# Eigenvectors 2 vs 3
axis1 <- 2
axis2 <- 3
plot(eigenvectors[,c(axis1,axis2)], xlab=paste0("Vector", axis1), ylab=paste0("Vector", axis2), type="n")
for(i in 1:length(populations)) {
  pop_i <- subset(allpops[,'sample'], allpops[,'super_pop'] %in% populations[i])
  eigen_i <- subset(eigenvectors, rownames(eigenvectors) %in% pop_i)
  points(eigen_i[,c(axis1,axis2)], col=colours[i])
}
points(eigen_outliers[,c(axis1,axis2)], col="turquoise", pch=19, cex=0.3)
#  legend("topleft", c(populations), col=colours, pch=21)

# Draw legend
plot(c(0,1),c(0,length(populations)),type='n',xlab='',ylab='')
x1<-rep(0,length(populations))
x2<-rep(0.2,length(populations))
y1=0:(length(populations)-1)
y2=1:(length(populations))
rect(x1,y1,x2,y2,col=colours)
for(i in 1:length(populations)) {
  text(0.4,(i-1)+0.5,populations[i])
}


write.table(ctrl_sd_outliers, paste0(out_prefix, "_outliers.txt"),quote=F,row.names=F,col.names=F)





# Write out the final set of Europeans:
fid <- NULL
iid <- NULL
for(i in 1:nrow(eigen_sample)) {
  iid <- c(iid, rownames(eigen_sample)[i])
  fid <- c(fid, subset(samples[,1], samples[,2] %in% rownames(eigen_sample)[i]))
}
finalEURset <- cbind(FID=fid, IID=iid)

write.table(finalEURset, paste0(out_prefix,"_finalEURset.txt"), quote=F, row.names=F, col.names=F)
