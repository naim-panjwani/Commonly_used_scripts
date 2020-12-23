setwd("/Users/naim/Desktop/1805-BIOJUME-Ethnicity/All_JME_and_1KG/")
# setwd("/Users/naim/Desktop/1805-BIOJUME-Ethnicity/All_JME_and_HGDP/")
# setwd("/Users/naim/Desktop/1805-BIOJUME-Ethnicity/All_JME_plus_OSC/")
# setwd("/Users/naim/Desktop/1805-BIOJUME-Ethnicity/All_JME_plus_1KG-EUR/")

# args=(commandArgs(TRUE))
# 
# eigenvalues_filename <- as.character(args[1]) # no header; two columns: 1st column is ignored and 2nd column has the eigenvectors
# eigenvectors_filename <- as.character(args[2]) # no header
# evec_iid_column <- as.integer(args[3])
# evec_1_column <- as.integer(args[4])
# evec_2_column <- as.integer(args[5])
# evec_3_column <- as.integer(args[6])
# iid_labels_file <- as.character(args[7]) # FID, IID, Label/ethnicity columns only; must have "sample" label; no header
# outprefix <- as.character(args[8]) # optional

eigenvalues_filename <- "34_eigenvalues.txt"
eigenvectors_filename <- "34_eigenvectors.txt"
evec_iid_column <- 1
evec_1_column <- 2
evec_2_column <- 3
evec_3_column <- 4
iid_labels_file <- "iid_labels.txt"
outprefix <- "35"

sphereSize = 12
sampleSphereSize = 16

if(is.na(outprefix) | is.null(outprefix)) outprefix <- ""

library(scales)
library(randomcoloR)
library(rgl)
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

M <- par3d("userMatrix")
st<-function(n,t)  {
  # to auto-rotate 3d plot 
  # eg. st(20,2)
  for (i in 1:n) {
    play3d( par3dinterp( userMatrix=list(M,rotate3d(M, pi, 1, 0, 0) ) ), duration=t )
    play3d( par3dinterp( userMatrix=list(M,rotate3d(M, pi, 0, 1, 1) ) ), duration=t)
    play3d( par3dinterp( userMatrix=list(M,rotate3d(M, pi, 0, 1, 0) ) ), duration=t)
    play3d( par3dinterp( userMatrix=list(M,rotate3d(M, pi, 1, 0, 1) ) ), duration=t)
  }
}
#====================================================================
# END OF FUNCTIONS
#====================================================================

# Plot screeplot and individuals' eigenvectors

# read files
eigenvalues <- read.table(eigenvalues_filename, header=F, stringsAsFactors = F)
eigenvectors <- read.table(eigenvectors_filename, header=F, stringsAsFactors = F)
iid_labels <- read.table(iid_labels_file, header=F, stringsAsFactors = F)
sample_individuals <- subset(iid_labels, iid_labels[,3] %in% c("sample","Sample"))
nonsample_individuals <- subset(iid_labels, !(iid_labels[,3] %in% c("sample","Sample")))


#### SCREEPLOT
#pdf(paste0(outprefix, "_screeplot.pdf"))
plot(eigenvalues[,2], xlim=c(1,nrow(eigenvalues)), ylab="Eigenvalue",xlab="Component number", main="Screeplot") #screeplot
#dev.off()


#### Draw box

# subset ethnic backgrounds and sample individuals
eigen_sample <- subset(eigenvectors, eigenvectors[,evec_iid_column] %in% sample_individuals[,2])
eigen_nonsample <- subset(eigenvectors, eigenvectors[,1] %in% nonsample_individuals[,2])
populations <- levels(as.factor(nonsample_individuals[,3]))
colours <- distinctColorPalette(length(populations)+1)

#draw the sample and nonsample individuals' first 3 PC points
plot3d(eigenvectors[,evec_1_column],
       eigenvectors[,evec_2_column],
       eigenvectors[,evec_3_column], type="n", xlab="E1", ylab="E2", zlab="E3")
for (i in 1:length(populations)){
  pop_i <- subset(nonsample_individuals, nonsample_individuals[,3] %in% populations[i])
  eigen_i <- subset(eigen_nonsample, eigen_nonsample[,evec_iid_column] %in% pop_i[,2])
points3d(eigen_i[,evec_1_column],
         eigen_i[,evec_2_column],
         eigen_i[,evec_3_column], col=colours[i], size=sphereSize)
}
points3d(eigen_sample[,evec_1_column],
         eigen_sample[,evec_2_column],
         eigen_sample[,evec_3_column],col=colours[length(populations)+1], size=sampleSphereSize)

# identify3d(eigenvectors[,evec_1_column],
#            eigenvectors[,evec_2_column],
#            eigenvectors[,evec_3_column],labels=eigenvectors[,evec_iid_column])   #### iid labels

# identify3d(eigenvectors[,evec_1_column],
#            eigenvectors[,evec_2_column],
#            eigenvectors[,evec_3_column],labels=iid_labels[match(eigenvectors[,evec_iid_column], iid_labels[,2]),3])  #### ethnic labels


# Draw the legend
n=length(populations)+1
par(mfrow=c(1,1))
plot(c(0,1),c(0,n),type='n',xlab='',ylab='')
x1<-rep(0,n)
x2<-rep(0.2,n)
y1=0:(n-1)
y2=1:n
rect(x1,y1,x2,y2,col=colours)
x=0.4
y=y2-0.5
pops <- c(populations, "Samples")
for(i in 1:n){
  text(x, y[i], pops[i])
}


# Want to determine ethnicity of each individual
pop_i_samples <- list()
esample_individuals <- eigen_sample
print(paste("Number of individuals per population"))
populations <- levels(as.factor(nonsample_individuals[,3]))
# populations[which(populations %in% "EUR")] <- populations[1] # want to classify the EUR samples before AMR as AMR overlaps EUR
# populations[1] <- "EUR"
for (i in 1:length(populations)){
  pop_i <- subset(nonsample_individuals, nonsample_individuals[,3] %in% populations[i])
#  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_nonsample, eigen_nonsample[,evec_iid_column] %in% pop_i[,2])
  individual_indices <- isInBox(eigen_i[,c(evec_1_column,evec_2_column,evec_3_column)],esample_individuals[,c(evec_1_column,evec_2_column,evec_3_column)]) 
  print(paste(populations[i],sum(individual_indices)))
  pop_i_samples[[i]] <- as.character(esample_individuals[individual_indices, evec_iid_column])
  names(pop_i_samples)[i] <- populations[i]
  if (sum(individual_indices)>0) esample_individuals <- esample_individuals[-which(individual_indices),] # remove sample individuals to prevent duplicates
}
sapply(pop_i_samples,length)
sum(sapply(pop_i_samples,length))

# who has not been classified?
esample_individuals[,evec_iid_column]
dim(esample_individuals)
with(esample_individuals, points3d(V2,V3,V4, col="grey", size=sampleSphereSize+3, pch="*"))

# include last category of unclassified individuals:
pop_i_samples[[length(pop_i_samples)+1]] <- as.character(esample_individuals[,evec_iid_column]) 
# populations <- c(populations, "Unknown")
names(pop_i_samples)[length(pop_i_samples)] <- "Unknown"


# Convert pop_i_samples to a data.frame:
ethnicity<-NULL
fid<-NULL
iid<-NULL
famfile <- iid_labels
getFID <- function(IID) {
  fid <- NULL
  for(i in 1:length(IID)) {
    fid <- c(fid, famfile[which(famfile[,2] %in% IID[i]),1])
  }
  return(fid)
}
for (i in 1:length(pop_i_samples)) {
  ethnicity <- c(ethnicity,rep(names(pop_i_samples)[i], length(pop_i_samples[[i]])))
  fid <- c(fid, getFID(pop_i_samples[[i]]))
  iid <- c(iid, pop_i_samples[[i]])
}
ind_ethnicities <- data.frame(FID=fid, IID=iid, Ethnicity=ethnicity)
#write.table(ind_ethnicities, paste0(outprefix,"ethnicities.txt"),
#            quote=F, row.names=F, col.names=T)

# Grab EUR, AMR and Unknown ethnic individuals for next PCA round



# getFID <- function(x) return(gsub("(.*):.*","\\1", x))
# getIID <- function(x) return(gsub(".*:(.*)","\\1", x))

