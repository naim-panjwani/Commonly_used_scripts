# This code is older than the "manual" ethnic analysis R code files
# This one was used for Eigenstrat output and is modified from practicum days
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
#====================================================================
# END OF FUNCTIONS
#====================================================================

# Plot screeplot and individuals' eigenvectors

# read files
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v3")

eigenvalues <- read.table("44_dataset_hapmap_autosomalSNPs_DI_rm.eigenvalue",header=FALSE)
eigenvectors <- read.table("44_dataset_hapmap_autosomalSNPs_DI_rm.eigenvector")
# outliers <- paste( read.table("28_ceu_outliers",header=F)$V1,
#                    read.table("28_ceu_outliers",header=F)$V2,
#                    sep=":")

# DO NOT ENTER THE MERGED DATASET; ONLY DATASET WITH NON-REFERENCE INDIVIDUALS
sample_individuals <- paste( read.table("35_dataset_flipped_snps.fam")$V1,
                             read.table("35_dataset_flipped_snps.fam")$V2,
                             sep=":")

hapmap_individuals <- paste(read.table("/Users/naim/Documents/Strug/hapmap_data/07_hapmap_hg19_related_rm_ceu_tsi.fam")$V1,
                            read.table("/Users/naim/Documents/Strug/hapmap_data/07_hapmap_hg19_related_rm_ceu_tsi.fam")$V2,
                            sep=":")
relationships <- read.table("hapmap3r3_b36_relationship.txt", header=TRUE)
relationships$ID <- paste(relationships$FID, relationships$IID, sep=":")

#### SCREEPLOT
#png("20_screeplot.png")
plot(eigenvalues$V1, xlim=c(1,10), ylab="Eigenvalue",xlab="Component number", main="Screeplot") #screeplot
#dev.off()


#### Draw box
library(rgl)

# subset ethnic backgrounds and sample individuals
eigen_sample <- subset(eigenvectors, eigenvectors$V1 %in% sample_individuals)
eigen_hapmap <- subset(eigenvectors, eigenvectors$V1 %in% relationships$ID)
#eigen_outliers <- subset(eigenvectors, eigenvectors$V1 %in% outliers)
populations <- levels(relationships$population)
colours <- c("red", "green", "purple", "black", "brown", "blue", 
             "cadetblue", "pink", 54, 53, 55, "orange")

#draw the sample and hapmap individuals' first 3 PC points
with(eigenvectors, plot3d(V2,V3,V4, type="n"))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_hapmap, eigen_hapmap$V1 %in% pop_i$ID)
  with(eigen_i, points3d(V2,V3,V4, col=colours[i+1]))
}
with(eigen_sample, points3d(V2,V3,V4,col=colours[1], size=6))

# Read in outliers:
outliers <- paste(read.table("46_outliers.txt",header=F)$V1,
                  read.table("46_outliers.txt",header=F)$V2,
                  sep=":")
outliers_removed <- subset(sample_individuals, !(sample_individuals %in% outliers))
#==========================================================================================================================
vectorIIDsTOlist <- function(samples) {
  fid <- iid <- NULL
  for (i in 1:length(samples)) {
    fid <- c(fid, gsub("(.*):.*","\\1",samples[i]))
    iid <- c(iid, gsub(".*:(.*)","\\1",samples[i]))
  }
  return(data.frame(FID=fid,IID=iid))
}
#==========================================================================================================================

# Draw the legend
par(mfrow=c(1,1))
plot(c(0,1),c(0,12),type='n',xlab='',ylab='')
x1<-rep(0,12)
x2<-rep(0.2,12)
y1=0:11
y2=1:12
rect(x1,y1,x2,y2,col=colours)
#x.text=rep(0.4,12)
#y.text=1:12
#text(x.text,y.text,express
text(0.4,0.5,expression('Samples'))
text(0.4,1.5,expression('ASW-African'))
text(0.4,2.5,expression('CEU-Caucasian'))
text(0.4,3.5,expression('CHB-Chinese-Beijing'))
text(0.4,4.5,expression('CHD-Chinese-Denver'))
text(0.4,5.5,expression('GIH-Indian'))
text(0.4,6.5,expression('JPT-Japanese'))
text(0.4,7.5,expression('LWK-Luhya'))
text(0.4,8.5,expression('MEX-Mexican'))
text(0.4,9.5,expression('MKK-Maasia'))
text(0.4,10.5,expression('TSI-Italian'))
text(0.4,11.5,expression('YRI-Yoruba'))


# Want to determine ethnicity of each individual
pop_i_samples <- list()
esample_individuals <- eigen_sample
print(paste("Number of individuals per population"))
populations <- levels(relationships$population)
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_hapmap, eigen_hapmap$V1 %in% pop_i$ID)
  individual_indices <- isInBox(eigen_i[,2:4],esample_individuals[,2:4]) 
  print(paste(populations[i],sum(individual_indices)))
  pop_i_samples[[i]] <- as.character(esample_individuals$V1[individual_indices])
  names(pop_i_samples)[i] <- populations[i]
  if (sum(individual_indices)>0) esample_individuals <- esample_individuals[-which(individual_indices),] # remove sample individuals to prevent duplicates
}
sapply(pop_i_samples,length)
sum(sapply(pop_i_samples,length))

# who has not been classified?
esample_individuals$V1
dim(esample_individuals)
with(esample_individuals, points3d(V2,V3,V4, col="grey", size=8, pch="*"))

# include last category of unclassified individuals:
pop_i_samples[[12]] <- as.character(esample_individuals$V1) 
populations <- c(populations, "Unknown")
names(pop_i_samples)[12] <- populations[12]


# Convert pop_i_samples to a data.frame:
ethnicity<-NULL
fid<-NULL
iid<-NULL
for (i in 1:length(pop_i_samples)) {
  ethnicity <- c(ethnicity,rep(names(pop_i_samples)[i], length(pop_i_samples[[i]])))
  fid <- c(fid, gsub("(.*):.*","\\1",pop_i_samples[[i]]))
  iid <- c(iid, gsub(".*:(.*)","\\1",pop_i_samples[[i]]))
}
ind_ethnicities <- data.frame(FID=fid, IID=iid, Ethnicity=ethnicity)
write.table(ind_ethnicities, "ethnicities.txt",
            quote=F, row.names=F, col.names=T)



getFID <- function(x) return(gsub("(.*):.*","\\1", x))
getIID <- function(x) return(gsub(".*:(.*)","\\1", x))



### List the sample individuals plus Hapmap individuals to run
### smartpca on to check for outliers more than 6 stdevs away

# Outliers away from the CEU/TSI cluster
# 1. Remove sample individuals that fall into other clusters:
ceu_tsi_sample_ind <- subset(ind_ethnicities, ind_ethnicities$Ethnicity %in% c("CEU", "TSI", "Unknown"))
dim(ceu_tsi_sample_ind)
# 134

# 2. List Hapmap CEU and TSI individuals only:
ceu_tsi_hapmap_ind <- subset(relationships, relationships$population %in% c("CEU","TSI"))[,c(8,7)]
hapmap_ceu_tsi_fids <- as.character(sapply(ceu_tsi_hapmap_ind$ID, getFID))
hapmap_ceu_tsi_iids <- as.character(sapply(ceu_tsi_hapmap_ind$ID, getIID))
ceu_tsi_hapmap_ind <- cbind(hapmap_ceu_tsi_fids, hapmap_ceu_tsi_iids, as.character(ceu_tsi_hapmap_ind[,2]))
colnames(ceu_tsi_hapmap_ind) <- c("FID","IID","Ethnicity")
dim(ceu_tsi_hapmap_ind)
# 267

# 3. Merge them for extraction using PLINK and subsequent PCA to detect individuals 6 stdevs away
ceu_merge <- rbind(ceu_tsi_sample_ind, ceu_tsi_hapmap_ind)
ceu_ind_for_pca <- write.table(ceu_merge, "ceu_ind_for_pca.txt",
                               quote=F, row.names=F, col.names=F)




# 4. After identifying the outliers, create final CEU/TSI list of individuals to keep and analyze
ceu_tsi_sample_ids <- paste(ceu_tsi_sample_ind$FID,ceu_tsi_sample_ind$IID,sep=":")
final_ceu <- subset(ceu_tsi_sample_ind, !(ceu_tsi_sample_ids %in% outliers[1:15]))
unknown_CEU_ind <- subset(final_ceu,final_ceu$Ethnicity %in% "Unknown")
# 26 Unknown individuals classified as CEU

final_ceu$Ethnicity <- as.character(final_ceu$Ethnicity)
final_ceu$Ethnicity[which(final_ceu$Ethnicity %in% "Unknown")] <- "CEU-TSI"

write.table(final_ceu, "ceu_individuals.txt",
            quote=F, row.names=F, col.names=F)

