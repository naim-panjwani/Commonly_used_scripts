
#setwd("/Users/naim/Documents/Strug/141202-Spit_for_Science_all_samples/QC")

eigenvalues <- read.table("36_eigenvalues.txt",header=FALSE)
eigenvectors <- read.table("35_mergedset_pcapc.ped")
case_individuals <- read.table("22_het_hap_snps_rm.fam")[,c(1:2)]
ctrl_individuals <- read.table("hapmap_data/07_hapmap_hg19_related_rm_ceu_tsi.fam")[,c(1:2)]

#### SCREEPLOT
png("36_screeplot.png")
plot(eigenvalues$V1, xlim=c(1,10), ylab="Eigenvalue",xlab="Component number", main="Screeplot") #screeplot
dev.off()


#library(rgl)


getFID <- function(IIDs) {
  fid <- character(length(IIDs))
  for (i in 1:length(IIDs)) {
    fid[i] <- as.character(eigenvectors$V1[which(as.character(eigenvectors$V2) %in% as.character(IIDs[i]))])
  }
  return(fid)
}


# subset ethnic backgrounds and sample individuals
eigen_cases <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(case_individuals$V2))
eigen_ctrls <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(ctrl_individuals$V2))

#draw the case and control individuals' first 3 PC points
# with(eigenvectors, plot3d(V7,V8,V9, type="n"))
# with(eigen_cases, points3d(V7,V8,V9,col="red", size=6))
# with(eigen_ctrls, points3d(V7,V8,V9,col="blue", size=6))


# 2D PLOTS
#png("19_pca_2d_plots.png")
#postscript("13_pca_2d_plots.ps", width=3, height=50)
# par(mfrow=c(3,1))
# # Vector 1 vs 2
# with(eigenvectors, plot(V7,V8,type="n", xlab="Vector1", ylab="Vector2"))
# with(eigen_cases, points(V7,V8, col="red"))
# with(eigen_ctrls, points(V7,V8, col="blue"))
# 
# 
# # Vector 2 vs 3
# with(eigenvectors, plot(V8,V9,type="n", xlab="Vector2", ylab="Vector3"))
# with(eigen_cases, points(V8,V9, col="red"))
# with(eigen_ctrls, points(V8,V9, col="blue"))
# 
# # Vector 1 vs 3
# with(eigenvectors, plot(V7,V9,type="n", xlab="Vector1", ylab="Vector3"))
# with(eigen_cases, points(V7,V9, col="red"))
# with(eigen_ctrls, points(V7,V9, col="blue"))
# #dev.off()
# par(mfrow=c(1,1))



# Identify outliers > 6 standard deviations away along any of the first four eigenvectors:
# (left_limit <- mean(eigenvectors$V7) - 6*sd(eigenvectors$V7))
# (right_limit <- mean(eigenvectors$V7) + 6*sd(eigenvectors$V7))
# (outliers_vec1 <- subset(eigenvectors, eigenvectors$V7 < left_limit | eigenvectors$V7 > right_limit))
# 
# (left_limit <- mean(eigenvectors$V8) - 6*sd(eigenvectors$V8))
# (right_limit <- mean(eigenvectors$V8) + 6*sd(eigenvectors$V8))
# (outliers_vect2 <- subset(eigenvectors, eigenvectors$V8 < left_limit | eigenvectors$V8 > right_limit))
# 
# (left_limit <- mean(eigenvectors$V9) - 6*sd(eigenvectors$V9))
# (right_limit <- mean(eigenvectors$V9) + 6*sd(eigenvectors$V9))
# (outliers_vect3 <- subset(eigenvectors, eigenvectors$V9 < left_limit | eigenvectors$V9 > right_limit))
# 
# (left_limit <- mean(eigenvectors$V10) - 6*sd(eigenvectors$V10))
# (right_limit <- mean(eigenvectors$V10) + 6*sd(eigenvectors$V10))
# (outliers_vect4 <- subset(eigenvectors, eigenvectors$V10 < left_limit | eigenvectors$V10 > right_limit))
# 
# outliers_iid <- unique(c(as.character(outliers_vec1$V2), as.character(outliers_vect2$V2), as.character(outliers_vect3$V2), as.character(outliers_vect4$V2)))
# outliers_fid <- getFID(outliers_iid)
# 
# (outliers <- data.frame(FID=outliers_fid, IID=outliers_iid))
#write.table(outliers, "13_pca_outliers.txt",quote=F,row.names=F,col.names=F)





# Identify outliers 6 std.devs away from Hapmap CEU/TSI clusters
# Identify outliers > 6 standard deviations away along any of the first four eigenvectors:
(left_limit <- mean(eigen_ctrls$V7) - 6*sd(eigen_ctrls$V7))
(right_limit <- mean(eigen_ctrls$V7) + 6*sd(eigen_ctrls$V7))
(outliers_vec1 <- subset(eigenvectors, eigenvectors$V7 < left_limit | eigenvectors$V7 > right_limit))

(left_limit <- mean(eigen_ctrls$V8) - 6*sd(eigen_ctrls$V8))
(right_limit <- mean(eigen_ctrls$V8) + 6*sd(eigen_ctrls$V8))
(outliers_vect2 <- subset(eigenvectors, eigenvectors$V8 < left_limit | eigenvectors$V8 > right_limit))

(left_limit <- mean(eigen_ctrls$V9) - 6*sd(eigen_ctrls$V9))
(right_limit <- mean(eigen_ctrls$V9) + 6*sd(eigen_ctrls$V9))
(outliers_vect3 <- subset(eigenvectors, eigenvectors$V9 < left_limit | eigenvectors$V9 > right_limit))

(left_limit <- mean(eigen_ctrls$V10) - 6*sd(eigen_ctrls$V10))
(right_limit <- mean(eigen_ctrls$V10) + 6*sd(eigen_ctrls$V10))
(outliers_vect4 <- subset(eigenvectors, eigenvectors$V10 < left_limit | eigenvectors$V10 > right_limit))


outliers_iid <- unique(c(as.character(outliers_vec1$V2), as.character(outliers_vect2$V2), as.character(outliers_vect3$V2), as.character(outliers_vect4$V2)))
outliers_fid <- getFID(outliers_iid)

(outliers <- data.frame(FID=outliers_fid, IID=outliers_iid))
write.table(outliers, "36_pca_outliers.txt",quote=F,row.names=F,col.names=F)



eigen_outliers <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(outliers$IID))
# 2D PLOTS
png("36_pca_2d_plots.png")
#postscript("33_pca_2d_plots.ps", width=3, height=50)
par(mfrow=c(3,1))
# Vector 1 vs 2
with(eigenvectors, plot(V7,V8,type="n", xlab="Vector1", ylab="Vector2"))
with(eigen_cases, points(V7,V8, col="red"))
with(eigen_ctrls, points(V7,V8, col="blue"))
with(eigen_outliers, points(V7,V8, col="green",pch="*"))

# Vector 2 vs 3
with(eigenvectors, plot(V8,V9,type="n", xlab="Vector2", ylab="Vector3"))
with(eigen_cases, points(V8,V9, col="red"))
with(eigen_ctrls, points(V8,V9, col="blue"))
with(eigen_outliers, points(V8,V9, col="green",pch="*"))

# Vector 1 vs 3
with(eigenvectors, plot(V7,V9,type="n", xlab="Vector1", ylab="Vector3"))
with(eigen_cases, points(V7,V9, col="red"))
with(eigen_ctrls, points(V7,V9, col="blue"))
with(eigen_outliers, points(V7,V9, col="green",pch="*"))
dev.off()
par(mfrow=c(1,1))



