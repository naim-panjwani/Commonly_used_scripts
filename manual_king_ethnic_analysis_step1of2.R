# R program to identify CEU/TSI individuals according to the following method/rule:
# 1. Get everyone that is within 6 sd's of Hapmap's CEU/TSI clusters
# 2. Assign the certain CEU/TSI individuals with the "box" method (do not expand the box by N*sd)
# 3. Among those sample individuals who are >6 sd's , assign them to other ethnicities 
#    with the "box" method (do not expand the box by N*sd) 
# 4. Those who are successfully identified as belonging to another ethnic group get 
#    removed from the consideration list
# 5. Further examine the remaining individuals (CEU/TSI 6 sd individuals + unassigned ones) and 
#    remove any obvious outliers that may remain
# 6. Now take the Hapmap CEU/TSI + sample CEU/TSI (certain and uncertain) + remaining unassigned individuals, 
#    redo PCA with just the CEU/TSI Hapmap samples (rather than including all other ethnic populations) 
#    and take out those samples >6 sd away. Only do 1 round.



#setwd("/Users/naim/Documents/Strug/141202-Spit_for_Science_all_samples/QC/MANUAL_ETHNICITY")

eigenvalues <- read.table("31_eigenvalues.txt",header=FALSE)
eigenvectors <- read.table("31_mergedset_pcapc.ped")
case_individuals <- read.table("../22_het_hap_snps_rm.fam")[,c(1:2)]
ctrl_individuals <- read.table("/Users/naim/Documents/Strug/hapmap_data/06_hapmap_hg19_related_rm_hh_rm.fam")[,c(1:2)]
relationships <- read.table("/Users/naim/Documents/Strug/hapmap_data/hapmap3r3_b36_relationship.txt",header=T)

#### SCREEPLOT
png("32_screeplot.png")
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
    left_limit <- mean(eigen_sd_set[,i]) - 6*sd(eigen_sd_set[,i])
    right_limit <- mean(eigen_sd_set[,i]) + 6*sd(eigen_sd_set[,i])
    outliers_vec_i <- subset(eigen_set, eigen_set[,i] < left_limit | eigen_set[,i] > right_limit)
    
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
eigen_cases <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(case_individuals$V2))
eigen_ctrls <- subset(eigenvectors, as.character(eigenvectors$V2) %in% as.character(ctrl_individuals$V2))

nrow(eigen_cases)
#4505

# 1. Get everyone that is within 6 sd's of Hapmap's CEU/TSI clusters (var eigen_ceu_6sd)
pop_i <- subset(relationships, relationships$population %in% "CEU" | relationships$population %in% "TSI")
eigen_ceu <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
ceu_outliers <- getOutliers(eigen_cases, eigen_ceu, fid_column=1, iid_column=2, eigenvalue_columns=7:9)
eigen_outliers <- subset(eigen_cases, as.character(eigen_cases$V2) %in% as.character(ceu_outliers$IID))
eigen_ceu_6sd <- subset(eigen_cases, !(as.character(eigen_cases$V2) %in% as.character(ceu_outliers$IID)))
nrow(eigen_ceu_6sd)
#4442

# 2. Assign the certain CEU/TSI individuals with the "box" method (do not expand the box by N*sd)
# var eigen_sample_ceu
eigen_sample_ceu <- subset(eigen_ceu_6sd, isInBox(eigen_ceu[,7:9], eigen_ceu_6sd[,7:9]))
nrow(eigen_sample_ceu)
# 3939 (out of 4442)
eigen_not_so_ceu <- subset(eigen_ceu_6sd, !(isInBox(eigen_ceu[,7:9], eigen_ceu_6sd[,7:9])))

# 3. Among those sample individuals who are >6 sd's (var eigen_outliers), assign them to other ethnicities 
#    with the "box" method (do not expand the box by N*sd) 
# var pop_i_samples
pop_i_samples <- list()
esample_individuals <- eigen_outliers
print(paste("Number of individuals per population"))
populations <- levels(relationships$population)
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  #pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
  individual_indices <- isInBox(eigen_i[,7:9],esample_individuals[,7:9]) 
  print(paste(populations[i],sum(individual_indices)))
  pop_i_samples[[i]] <- as.character(esample_individuals$V2[individual_indices])
  names(pop_i_samples)[i] <- populations[i]
  if (sum(individual_indices)>0) esample_individuals <- esample_individuals[-which(individual_indices),] # remove sample individuals to prevent duplicates
}
# put the CEU/TSI's in CEU group
pop_i_samples[[2]] <- as.character(eigen_sample_ceu[,2])
sapply(pop_i_samples,length)
# ASW  CEU  CHB  CHD  GIH  JPT  LWK  MEX  MKK  TSI  YRI 
#   0 3939    1    0    3    0    0    5    0    0    0 

sum(sapply(pop_i_samples,length))
# 3948 (out of 4442)


# 4. Those who are successfully identified as belonging to another ethnic group get 
#    removed from the consideration list
# remaining unassigned outliers: var esample_individuals (54 individuals)


# 5. Further examine the remaining individuals (CEU/TSI 6 sd individuals + unassigned ones) and 
#    remove any obvious outliers that may remain (ie. Say > 12? sd away)
# Graphing needed

populations <- levels(relationships$population)
colours <- c("red", "green", "purple", "black", "brown", "blue", 
             "cadetblue", "pink", 54, 53, 55, "orange")

# draw the case and control individuals' first 3 PC points
with(eigenvectors, plot3d(V7,V8,V9, type="n"))
with(eigen_cases, points3d(V7,V8,V9,col=colours[1]))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  #pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
  with(eigen_i, points3d(V7,V8,V9, col=colours[i+1]))
}


# Mark unassigned outlier individuals with "green" ring
with(esample_individuals, points3d(V7,V8,V9, col="green", size=6))

# add a colored "ring" on top of sample points for those who got assigned to an ethnic cluster
for(i in 1:length(populations)) {
  sample_i <- subset(eigen_cases, as.character(eigen_cases[,2]) %in% as.character(pop_i_samples[[i]]))
  with(sample_i, points3d(V7,V8,V9,col=colours[i+1],size=10))
}


# Obvious outliers (var eyeballing_outliers):
# Var with eyeballing_outliers (out of the unassigned outliers >6 sd's from Hapmap CEU/TSI) removed: new_esample_ind
# Eigenvector3 (V9) < 
# eyeballing_outliers <- subset(esample_individuals, esample_individuals$V9 < -0.04)
# with(eyeballing_outliers, points3d(V7,V8,V9, col="blue", size=8))
# new_esample_ind <- esample_individuals[-which(as.character(esample_individuals[,2]) %in% as.character(eyeballing_outliers[,2])),]
# ?? individuals
# Give them a blue ring on top of the green ring


# Eigenvector2 (v8) < -0.01
#eyeballing_outliers <- rbind(eyeballing_outliers, subset(new_esample_ind, new_esample_ind$V8 < -0.1))
eyeballing_outliers <- subset(esample_individuals, esample_individuals$V8 < -0.01)
with(eyeballing_outliers, points3d(V7,V8,V9, col="blue", size=8))
new_esample_ind <- esample_individuals[-which(as.character(esample_individuals[,2]) %in% as.character(eyeballing_outliers[,2])),]


# Eigenvector1 (v7) > 0.02
eyeballing_outliers <- rbind(eyeballing_outliers, subset(new_esample_ind, new_esample_ind$V7 > 0.02))
with(eyeballing_outliers, points3d(V7,V8,V9, col="blue", size=8))
new_esample_ind <- esample_individuals[-which(as.character(esample_individuals[,2]) %in% as.character(eyeballing_outliers[,2])),]


nrow(eyeballing_outliers)
# 3
nrow(new_esample_ind)
# 51


# 51 individuals left who are of mixed ethnic background and, by eye, uncertain what to do


# Assign and save the ethnicities for samples that are certain
# include last category of unclassified individuals:
pop_i_samples[[12]] <- c(as.character(esample_individuals$V2), as.character(eigen_not_so_ceu[,2]))
populations <- c(populations, "Unknown")
names(pop_i_samples)[12] <- populations[12]
sapply(pop_i_samples,length)
# ASW     CEU     CHB     CHD     GIH     JPT     LWK     MEX     MKK     TSI     YRI Unknown 
#   0    3939       1       0       3       0       0       5       0       0       0     557 
sum(sapply(pop_i_samples,length))
# 4505



# Convert pop_i_samples to a data.frame and save:
ethnicity<-NULL
fid<-NULL
iid<-NULL
for (i in 1:length(pop_i_samples)) {
  ethnicity <- c(ethnicity,rep(names(pop_i_samples)[i], length(pop_i_samples[[i]])))
  fid <- c(fid, getFID(pop_i_samples[[i]]))
  iid <- c(iid, pop_i_samples[[i]])
}
ind_ethnicities <- data.frame(FID=fid, IID=iid, Ethnicity=ethnicity)
write.table(ind_ethnicities, "32_ethnicities.txt",
            quote=F, row.names=F, col.names=T)



# 2D plots
# Vector 1 vs 2
png("32_2d_pca_plots.png", width=1440,height=1440, pointsize=24)
par(mfrow=c(2,2))
with(eigenvectors, plot(V7,V8,type="n", xlab="Vector1", ylab="Vector2"))
with(eigen_cases, points(V7,V8, col=colours[1]))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  eigen_i <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
  with(eigen_i, points(V7,V8, col=colours[i+1]))
}
# Mark unassigned outlier individuals with "green" ring
with(esample_individuals, points(V7,V8, col="green", cex=1.5))
# add a colored "ring" on top of sample points for those who got assigned to an ethnic cluster
for(i in 1:length(populations)) {
  sample_i <- subset(eigen_cases, as.character(eigen_cases[,2]) %in% as.character(pop_i_samples[[i]]))
  with(sample_i, points(V7,V8,col=colours[i+1],cex=2.5))
}
# add a blue ring to the outliers you are "eyeballing"
with(eyeballing_outliers, points(V7,V8, col="blue", cex=2))

# Vector 2 vs 3
with(eigenvectors, plot(V8,V9,type="n", xlab="Vector2", ylab="Vector3"))
with(eigen_cases, points(V8,V9, col=colours[1]))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  eigen_i <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
  with(eigen_i, points(V8,V9, col=colours[i+1]))
}
# Mark unassigned outlier individuals with "green" ring
with(esample_individuals, points(V8,V9, col="green", cex=1.5))
# add a colored "ring" on top of sample points for those who got assigned to an ethnic cluster
for(i in 1:length(populations)) {
  sample_i <- subset(eigen_cases, as.character(eigen_cases[,2]) %in% as.character(pop_i_samples[[i]]))
  with(sample_i, points(V8,V9,col=colours[i+1],cex=2.5))
}
# add a blue ring to the outliers you are "eyeballing"
with(eyeballing_outliers, points(V8,V9, col="blue", cex=2))

# Vector 1 vs 3
with(eigenvectors, plot(V7,V9,type="n", xlab="Vector1", ylab="Vector3"))
with(eigen_cases, points(V7,V9, col=colours[1]))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  eigen_i <- subset(eigen_ctrls, as.character(eigen_ctrls$V2) %in% as.character(pop_i$IID))
  with(eigen_i, points(V7,V9, col=colours[i+1]))
}
# Mark unassigned outlier individuals with "green" ring
with(esample_individuals, points(V7,V9, col="green", cex=1.5))
# add a colored "ring" on top of sample points for those who got assigned to an ethnic cluster
for(i in 1:length(populations)) {
  sample_i <- subset(eigen_cases, as.character(eigen_cases[,2]) %in% as.character(pop_i_samples[[i]]))
  with(sample_i, points(V7,V9,col=colours[i+1],cex=2.5))
}
# add a blue ring to the outliers you are "eyeballing"
with(eyeballing_outliers, points(V7,V9, col="blue", cex=2))



#draw legend
plot(c(0,1),c(0,16),type='n',xlab='',ylab='')
x1<-rep(0,12)
x2<-rep(0.2,12)
y1=4:15
y2=5:16
rect(x1,y1,x2,y2,col=colours)
#x.text=rep(0.2,12)
#y.text=1:12
#text(x.text,y.text,express
text(0.2,4.5,expression('Samples'),pos=4)
text(0.2,5.5,expression('ASW-African'),pos=4)
text(0.2,6.5,expression('CEU-Caucasian'),pos=4)
text(0.2,7.5,expression('CHB-Chinese-Beijing'),pos=4)
text(0.2,8.5,expression('CHD-Chinese-Denver'),pos=4)
text(0.2,9.5,expression('GIH-Indian'),pos=4)
text(0.2,10.5,expression('JPT-Japanese'),pos=4)
text(0.2,11.5,expression('LWK-Luhya'),pos=4)
text(0.2,12.5,expression('MEX-Mexican'),pos=4)
text(0.2,13.5,expression('MKK-Maasia'),pos=4)
text(0.2,14.5,expression('TSI-Italian'),pos=4)
text(0.2,15.5,expression('YRI-Yoruba'),pos=4)

points(0.1,0.5,cex=2,col="blue",pch=19)
points(0.1,0.5,cex=1.5,col="green",pch=19)
points(0.1,0.5,col="red",cex=1,pch=19)
text(0.2,0.5,expression('Eyeballed outlier'),pos=4)

points(0.9,0.5,cex=2,col="blue")
points(0.9,0.5,cex=1.5,col="green")
points(0.9,0.5,col="red",cex=1)

points(0.1,1.5,cex=1.5,col="green",pch=19)
points(0.1,1.5,col="red",cex=1,pch=19)
text(0.2,1.5,expression('>6 SD away from Hapmap CEU/TSI'),pos=4)

points(0.9,1.5,cex=1.5,col="green")
points(0.9,1.5,col="red",cex=1)

points(0.1,2.5,cex=2.5,col="pink",pch=19)
points(0.1,2.5,col="red",cex=1,pch=19)
text(0.2,2.5,expression('Assigned to the ethnic group'),pos=4)

points(0.9,2.5,cex=2.5,col="pink")
points(0.9,2.5,col="red",cex=1)

dev.off()



# 6. Now take the Hapmap CEU/TSI + sample CEU/TSI (certain and uncertain) + remaining unassigned individuals, 
#    redo PCA with just the CEU/TSI Hapmap samples (rather than including all other ethnic populations) 
#    and take out those samples >6 sd away. Only do 1 round.

# Hapmap CEU/TSI:
eigen_ceu
# Sample CEU/TSI + not-so-CEU:
eigen_ceu_6sd
# Remaining uknown individuals:
new_esample_ind

redoPCA_ind <- rbind(eigen_ceu, eigen_ceu_6sd, new_esample_ind)
# 4442+51=4493 sample individuals are going forward and 12 are essentially out due to assignment to another ethnicity or 
# because they have been eyeballed as being outliers

write.table(redoPCA_ind[,1:2], "32_redoPCAind.txt", quote=F, row.names=F, col.names=F)
