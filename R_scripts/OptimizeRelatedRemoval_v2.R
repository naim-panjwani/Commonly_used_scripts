# Reads in .genome and .imiss files, then 
# 1. generates siblings.txt, parent_offspring.txt, duplicates.txt, second_degree_relatives.txt and unrelated.txt (PI_HAT>=0.05) files
# 2. optimizes the list for removal of related individuals according to # of relationships and missingness


#library(rgl)
library(lattice)
#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/150304-GWAS_GEE_Genotyped_CTS_or_SSD_CORRECTION/simple_logistic_iid/")

# Read .genome and .imiss files
#relatedness <- read.table("17_Step12_kinship_IID_updated.genome", header=TRUE)
imiss <- read.table("CTS_or_SSD_CEU_dbSNP141.imiss", header=TRUE)
#with(relatedness, plot3d(Z0,Z1,Z2))

# relatedness<-subset(relatedness, relatedness$PI_HAT!=0)
# print(densityplot(~relatedness$PI_HAT, width=0.04, xlab="PI HAT")) 
# print(densityplot(~relatedness$PI_HAT, from=0.0, to=0.7,width=0.04, xlab="PI HAT"))
# print(densityplot(~relatedness$PI_HAT, from=0.0, to=0.3,width=0.04, xlab="PI HAT")) # some 3rd-degree relatives (ie. first cousins) but hard to resolve


# unrelated2 <- subset(relatedness, with(relatedness, Z2<0.23 & PI_HAT < 0.24))
#with(relatedness, plot3d(Z0,Z1,Z2, type="n"))
#with(unrelated2, points3d(Z0,Z1,Z2, col="black"))
# unrelated2_pairs <- paste(as.character(unrelated2$IID1), as.character(unrelated2$IID2),sep=",")
# first cousins considered unrelated

# second_degree_relatives <- subset(relatedness, with(relatedness, Z0>0.3 & Z2<0.1 & 
#                                                       PI_HAT >=0.24 & PI_HAT <0.4))
# second_degree_pairs <- paste(as.character(second_degree_relatives$IID1), as.character(second_degree_relatives$IID2),sep=",")
# #with(second_degree_relatives, points3d(Z0,Z1,Z2, col="darkgoldenrod3"))
# sib_pairs <- subset(relatedness, with(relatedness, Z0>0.1 & Z1>0.3 & Z2>=0.1))
# sib_pairs_pairs <- paste(as.character(sib_pairs$IID1), as.character(sib_pairs$IID2),sep=",")
# #with(sib_pairs, points3d(Z0,Z1,Z2, col="blue"))
# parent_offspring <- subset(relatedness, with(relatedness, Z0<0.1 & Z1>0.8 & Z2<0.2))
# parent_offspring_pairs <- paste(as.character(parent_offspring$IID1), as.character(parent_offspring$IID2),sep=",")
# #with(parent_offspring, points3d(Z0,Z1,Z2, col="red"))
# twins_duplicates <- subset(relatedness, with(relatedness, Z2>0.9))
# duplicate_pairs <- paste(as.character(twins_duplicates$IID1), as.character(twins_duplicates$IID2),sep=",")
# #with(twins_duplicates, points3d(Z0,Z1,Z2, col="darkgreen"))
# 
# (all_included <- dim(relatedness)[1] == (dim(unrelated2)[1]+dim(second_degree_relatives)[1]+dim(sib_pairs)[1]+
#                                            dim(parent_offspring)[1]+dim(twins_duplicates)[1]) )
# 
# write.table(second_degree_relatives, "18_second_degree_relatives.txt",
#             quote=FALSE,row.names=FALSE, col.names=TRUE)
# write.table(sib_pairs, "18_siblings.txt",
#             quote=FALSE,row.names=FALSE, col.names=TRUE)
# write.table(parent_offspring, "18_parent_offspring.txt",
#             quote=FALSE,row.names=FALSE, col.names=TRUE)
# write.table(twins_duplicates, "18_duplicates.txt",
#             quote=FALSE,row.names=FALSE, col.names=TRUE)
# write.table(unrelated2, "18_unrelated.txt",
#             quote=FALSE,row.names=FALSE, col.names=TRUE)
# 
# 
# 

WorstPair <- function(related_pairs, imiss) {
  # PRE: takes in .genome format data.frame and .imiss missingness information
  # POST: returns vector of each pair in .genome file with highest missingness
  
  num_pairs <- dim(related_pairs)[1]
  to_be_removed <- NULL
  for (i in 1:num_pairs) { # for each pair
    miss1 <- subset(imiss, imiss$IID %in% related_pairs$IID1[i])$F_MISS
    miss2 <- subset(imiss, imiss$IID %in% related_pairs$IID2[i])$F_MISS
    if(miss1 > miss2) { 
      pair_to_be_removed <- cbind(as.character(related_pairs$FID1[i]),as.character(related_pairs$IID1[i])) 
    } else {
      pair_to_be_removed <- cbind(as.character(related_pairs$FID2[i]),as.character(related_pairs$IID2[i])) 
    }
    to_be_removed <- rbind(to_be_removed,pair_to_be_removed)
  }
  return(to_be_removed)
}


countDuplicates <- function(x, v) {
  k<-0
  for (i in 1:length(v)) {
    if(x==v[i]) {
      k<-k+1
    }
  }
  return(k)
}
#------------------------------------------------------------------------------------------------
OptimizeRelatedRemoval <- function(related_pairs, imiss) {
  # Version 2
  # Initialize variables
  to_be_removed <- NULL
  fids <- NULL
  j <- 1
  
  related_pairs <- trim(related_pairs,imiss)
  
  # Identify unique IID vector
  individuals <- unique(c(as.character(related_pairs$IID1), as.character(related_pairs$IID2)))
  as_is_list <- c(as.character(related_pairs$IID1), as.character(related_pairs$IID2))
  count <- integer(length(individuals))  
  
  # For each individual, count the number of relationships
  for(i in 1:length(individuals)) { # for each individual
    count[i] <- length(which(as_is_list %in% individuals[i]))
  }
  
  for (i in 1:nrow(related_pairs)) { # for each pair
    iid1 <- as.character(related_pairs$IID1[i])
    iid2 <- as.character(related_pairs$IID2[i])
    if(!(iid1 %in% to_be_removed) & !(iid2 %in% to_be_removed)) {
      IID1_rel_count <- count[which(individuals %in% iid1)]
      IID2_rel_count <- count[which(individuals %in% iid2)]
      if(IID1_rel_count > IID2_rel_count) {
        to_be_removed[j] <- iid1
        fids[j]<-as.character(related_pairs$FID1[i])
        j<-j+1
      } else if (IID1_rel_count == IID2_rel_count) {
        temp <- WorstPair(related_pairs[i,],imiss)
        to_be_removed[j] <- temp[1,2]
        fids[j] <- temp[1,1]
        j<-j+1
      } else {
        to_be_removed[j] <- iid2
        fids[j]<-as.character(related_pairs$FID2[i])
        j<-j+1
      }
    }
  }
  return(as.data.frame(cbind(FID=as.character(fids),IID=as.character(to_be_removed))))
}
#------------------------------------------------------------------------------------------------
trim <- function(related, imiss) {
  # Removes individuals in "related" that are not present in "imiss"
  
  trimed_related <- NULL
  for(i in 1:nrow(related)) {
    iid1 <- as.character(related$IID1[i])
    iid2 <- as.character(related$IID2[i])
    available_set <- as.character(imiss$IID)
    if( (iid1 %in% available_set) & (iid2 %in% available_set) ) {
      trimed_related <- rbind(trimed_related, related[i,])
    }
  }
  return(trimed_related)
}
#------------------------------------------------------------------------------------------------

related <- rbind(twins_duplicates,parent_offsprings,sib_pairs,second_degree_pairs)
related_individuals_to_remove <- OptimizeRelatedRemoval(related, imiss)

write.table(related_individuals_to_remove, "related_individuals_to_remove.txt",
            quote=F, row.names=F, col.names=F)

