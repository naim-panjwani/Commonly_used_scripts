# R program to identify outliers >6SD's away and remove them:


#setwd("/Users/naim/Documents/Strug/141202-Spit_for_Science_all_samples/QC/MANUAL_ETHNICITY")

eigenvalues <- read.table("37_eigenvalues.txt",header=FALSE)
eigenvectors <- read.table("37_mergedset_pcapc.ped")
case_individuals <- read.table("../22_het_hap_snps_rm.fam")[,c(1:2)]
ctrl_individuals <- read.table("/Users/naim/Documents/Strug/hapmap_data/06_hapmap_hg19_related_rm_hh_rm.fam")[,c(1:2)]
relationships <- read.table("/Users/naim/Documents/Strug/hapmap_data/hapmap3r3_b36_relationship.txt",header=T)

#### SCREEPLOT
png("38_screeplot.png")
plot(eigenvalues$V1, xlim=c(1,10), ylab="Eigenvalue",xlab="Component number", main="Screeplot") #screeplot
dev.off()


library(rgl)

#==============================================================================================================
getFID <- function(IIDs) {
  fid <- character(length(IIDs))
  for (i in 1:length(IIDs)) {
    fid[i] <- as.character(eigenvectors$V1[which(as.character(eigenvectors$V2) %in% as.character(IIDs[i]))])
  }
  return(fid)
}
#-------------------
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
#------------------
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
#------------------
isInRange <- function(range, point) (point >= range[1] & point <= range[2])
#------------------
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
#------------------
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
#==============================================================================================================

# subset ethnic backgrounds and sample individuals
eigenvectors <- subset(eigenvectors, !(eigenvectors[,7]=="X"))
for(i in 7:ncol(eigenvectors)) eigenvectors[,i] <- as.numeric(as.character(eigenvectors[,i]))
eigen_cases <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(case_individuals$V2))
eigen_ctrls <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(ctrl_individuals$V2))

# Verify all controls are CEU/TSI
pop_i <- subset(relationships, relationships$population %in% "CEU" | relationships$population %in% "TSI")
eigen_hapmapceu <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
(nrow(eigen_ctrls)==nrow(eigen_hapmapceu))

# Outliers from the overall set:
outliers <- getOutliers(eigen_set=eigenvectors, eigen_sd_set=eigenvectors, fid_column=1, iid_column=2, eigenvalue_columns=7:9)
eigen_outliers <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(outliers$IID))
eigen_ceu_6sd <- subset(eigenvectors, !(as.character(eigenvectors$V2) %in% as.character(outliers$IID)))

nrow(eigen_ceu_6sd)
# 4706
nrow(eigen_outliers)
# 0

colours <- c("red", "purple")

# draw the case and control individuals' first 3 PC points
with(eigenvectors, plot3d(V7,V8,V9, type="n"))
with(eigen_cases, points3d(V7,V8,V9,col=colours[1]))
with(eigen_ctrls, points3d(V7,V8,V9,col=colours[2]))
# Mark outliers individuals with "green" ring
with(eigen_outliers, points3d(V7,V8,V9, col="green", size=8))


# 2D plots
# Vector 1 vs 2
png("38_2d_pca_plots.png", width=1440,height=1440, pointsize=24)
par(mfrow=c(2,2))
with(eigenvectors, plot(V7,V8,type="n", xlab="Vector1", ylab="Vector2"))
with(eigen_cases, points(V7,V8, col=colours[1]))
with(eigen_ctrls, points(V7,V8, col=colours[2]))
# Mark outlier individuals with "green" ring
with(eigen_outliers, points(V7,V8, col="green", cex=1.5))

# Vector 2 vs 3
with(eigenvectors, plot(V8,V9,type="n", xlab="Vector2", ylab="Vector3"))
with(eigen_cases, points(V8,V9, col=colours[1]))
with(eigen_ctrls, points(V8,V9, col=colours[2]))
# Mark outlier individuals with "green" ring
with(eigen_outliers, points(V8,V9, col="green", cex=1.5))


# Vector 1 vs 3
with(eigenvectors, plot(V7,V9,type="n", xlab="Vector1", ylab="Vector3"))
with(eigen_cases, points(V7,V9, col=colours[1]))
with(eigen_ctrls, points(V7,V9, col=colours[2]))
# Mark outlier individuals with "green" ring
with(eigen_outliers, points(V7,V9, col="green", cex=1.5))


#draw legend
plot(c(0,0.5),c(0,2),type='n',xlab='',ylab='')

points(0.1,0.5,cex=2,col="green")
points(0.1,0.5,col=colours[1],cex=1.5)
text(0.11,0.5,expression('Outliers'),pos=4)

points(0.1,1,col=colours[1],cex=1.5)
text(0.11,1,expression('Sample'),pos=4)

points(0.1,1.5,col=colours[2],cex=1.5)
text(0.11,1.5,expression('Hapmap CEU/TSI'),pos=4)

dev.off()


#final_set <- eigen_cases[-which(as.character(eigen_cases[,2]) %in% as.character(outliers[,2])),1:2]
final_set <- eigen_cases

write.table(final_set, "38_ceu_tsi_finalset.txt", quote=F, row.names=F, col.names=F)
