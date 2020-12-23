#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5")

parent_offsprings <- read.table("16_parent_offspring.txt",header=TRUE, stringsAsFactors = F)
sib_pairs <- read.table("16_siblings.txt", header=TRUE, stringsAsFactors = F)
second_degree_pairs <- read.table("16_second_degree_relatives.txt",header=TRUE, stringsAsFactors = F)
twins <- read.table("16_duplicates.txt",header=TRUE, stringsAsFactors = F)
unrelated2 <- read.table("16_unrelated.txt",header=TRUE, stringsAsFactors = F)
all_pairs <- read.table("16_Step12_kinship.genome",header=TRUE, stringsAsFactors = F)

fam <- read.table("Step1_5ind_rm_UKRE_sampleID.fam",header=FALSE, stringsAsFactors = F)

# Update IIDs First!

updateIIDs <- function(genome_file) {
  genome_file[,2] <- as.character(genome_file[,2])
  genome_file[,4] <- as.character(genome_file[,4])
  genome_file[which(as.character(genome_file[,2]) %in% "1035_2-1"),2] <- "1035_201"
  genome_file[which(as.character(genome_file[,2]) %in% "1059_301dup"),2] <- "1059_302"
  genome_file[which(as.character(genome_file[,2]) %in% "7015_301"),2] <- "7015_201"
  genome_file[which(as.character(genome_file[,2]) %in% "7015_301dup"),2] <- "7015_301"
  genome_file[which(as.character(genome_file[,2]) %in% "1061_302"),2] <- "1061_301"
  genome_file[which(as.character(genome_file[,4]) %in% "1035_2-1"),4] <- "1035_201"
  genome_file[which(as.character(genome_file[,4]) %in% "1059_301dup"),4] <- "1059_302"
  genome_file[which(as.character(genome_file[,4]) %in% "7015_301"),4] <- "7015_201"
  genome_file[which(as.character(genome_file[,4]) %in% "7015_301dup"),4] <- "7015_301"
  genome_file[which(as.character(genome_file[,4]) %in% "1061_302"),4] <- "1061_301"
  return(genome_file)
}

parent_offsprings <- updateIIDs(parent_offsprings)
sib_pairs <- updateIIDs(sib_pairs)
second_degree_pairs <- updateIIDs(second_degree_pairs)
twins <- updateIIDs(twins)
unrelated2 <- updateIIDs(unrelated2)
all_pairs <- updateIIDs(all_pairs)
#fam <- updateIIDs(fam)

second_degree_relatives<-second_degree_pairs
twins_duplicates <- twins

write.table(parent_offsprings, "17_parent_offspring_IID_updated.txt", quote=F, row.names=F,col.names=T)
write.table(sib_pairs, "17_siblings.txt_IID_updated.txt", quote=F, row.names=F,col.names=T)
write.table(second_degree_pairs, "17_second_degree_relatives_IID_updated.txt", quote=F, row.names=F,col.names=T)
write.table(twins, "17_duplicates_IID_updated.txt", quote=F, row.names=F,col.names=T)
write.table(unrelated2, "17_unrelated_IID_updated.txt", quote=F, row.names=F,col.names=T)
write.table(all_pairs, "17_Step12_kinship_IID_updated.genome", quote=F, row.names=F,col.names=T)


############################### FUNCTIONS #################################
pairs <- function(individual, related_pairs) {
  # PRE: individual is a character type and related_pairs has all the pairs 
  #      in two columns where individual should be present in at least one column
  # POST: returns all the individuals paired with individual in a character vector
  
  relatives <- NULL
  if (class(individual) != "character") stop("Individual ID must be of character type")
  else {
    for (i in 1:dim(related_pairs)[1]) {
      if(related_pairs[i,1]==individual | related_pairs[i,2]==individual) {
        ifelse(related_pairs[i,1]==individual,relatives<-c(relatives,as.character(related_pairs[i,2])),
               relatives<-c(relatives,as.character(related_pairs[i,1])))
      }
    }
  }
  return(relatives)
} # end of pairs
#----------------------------------------------------------------------------------------------
getFID <- function(IIDs, fam) {
  # PRE: IIDs is a vector of IID's contained within the fam file
  # POST: returns vector of FID's that matches the order of inputted IID's
  
  fids <- character(length(IIDs))
  for (i in 1:length(IIDs)) { # for each IID
    fids[i] <- as.character(fam[which(as.character(fam[,2])==as.character(IIDs[i])),1])
  }
  return(fids)
} # end of getFID
#----------------------------------------------------------------------------------------------
getSex <- function(IIDs, fam) {
  # PRE: IIDs is a vector of IID's contained within the fam file
  # POST: returns vector of IID's sex that matches the order of inputted IID's
  
  sex <- integer(length(IIDs))
  for (i in 1:length(IIDs)) { # for each IID
    sex[i] <- fam[which(as.character(fam[,2])==as.character(IIDs[i])),5]
  }
  return(sex)
} # end of getSex
#----------------------------------------------------------------------------------------------
permutePairs <- function(individuals){ # private function
  # PRE: give a vector of IDs
  # POST: permutes and returns all possible unique combinations
  
  n_perm <- sum(c((length(individuals)-1):1))
  pairs <- matrix(ncol=2,nrow=n_perm)
  if(length(individuals)==1) stop("Only one individual to permute")
  r<-1
  for (i in 1:(length(individuals)-1)) {
    for (j in (i+1):length(individuals)) {
      pairs[r,] <- c(individuals[i],individuals[j])
      r<-r+1
    }
  }
  return(pairs)
}
#----------------------------------------------------------------------------------------------
related <- function(individuals, rel_pairs){
  # PRE: individuals contains IIDs to be tested for relatedness as specified in rel_pairs;
  #      rel_pairs is two columns of IIDs of pairs of relatives
  # POST: returns TRUE if any of the given individuals are related as specified in rel_pairs
  
  rel <- FALSE
  possible_pairs <- permutePairs(individuals)
  for(i in 1:dim(possible_pairs)[1]) { # for each pair
    index<-which(as.character(rel_pairs[,1]) %in% possible_pairs[i,1])
    index2<-which(as.character(rel_pairs[,2]) %in% possible_pairs[i,1])
    if(length(index)>0) {
      for(j in 1:length(index)) {
        if(rel_pairs[index[j],2]==possible_pairs[i,2]) return(rel<-TRUE)
      }
    }
    if(length(index2)>0) {
      for(j in 1:length(index2)) {
        if(rel_pairs[index2[j],1]==possible_pairs[i,2]) return(rel<-TRUE)
      }
    }    
  }
  return(rel)
}
#----------------------------------------------------------------------------------------------
allRelated <- function(individuals, rel_pairs){
  # PRE: individuals contains IIDs to be tested for relatedness as specified in rel_pairs;
  #      rel_pairs is two columns of IIDs of pairs of relatives
  # POST: returns TRUE if ALL of the given individuals are related as specified in rel_pairs
  
  possible_pairs <- permutePairs(individuals)
  rel <- logical(dim(possible_pairs)[1])
  for(i in 1:dim(possible_pairs)[1]) { # for each pair
    index<-which(as.character(rel_pairs[,1]) %in% possible_pairs[i,1])
    index2<-which(as.character(rel_pairs[,2]) %in% possible_pairs[i,1])
    if(length(index)>0) {
      for(j in 1:length(index)) {
        if(rel_pairs[index[j],2]==possible_pairs[i,2]) rel[i]<-TRUE 
      }
    }
    if(length(index2)>0) {
      for(j in 1:length(index2)) {
        if(rel_pairs[index2[j],1]==possible_pairs[i,2]) rel[i]<-TRUE
      }
    }    
  }
  return(all(rel))
}
#----------------------------------------------------------------------------------------------
inRange <- function(value, Range) {
  return(ifelse(value>=Range[1] & value<=Range[2],TRUE,FALSE))
}
#----------------------------------------------------------------------------------------------
getRelationship_type <- function(z_values,pi_hat) {
  # PRE:
  # POST:
  
  result<-NULL
  if(inRange(pi_hat,range(sib_pairs$PI_HAT)) & inRange(z_values[2],range(sib_pairs$Z1))) {
    result<-"sib_pair"
  } else if(inRange(pi_hat,range(parent_offsprings$PI_HAT)) & inRange(z_values[2],range(parent_offsprings$Z1))) {
    result<-"parent_child"
  } else if(inRange(pi_hat,range(second_degree_relatives$PI_HAT))) {
    result<-"second_degree"
  } else if(inRange(pi_hat,range(twins_duplicates$PI_HAT))) {
    result<-"twin_duplicate"
  } else {
    result<-"unrelated"
  }
  return(result)
}
#----------------------------------------------------------------------------------------------
relationship_summary <- function(individuals, rel_pairs_detailed){
  # PRE: individuals contains IIDs to be tested for relatedness as specified in rel_pairs_detailed;
  #      rel_pairs is the genome file PLINK output with IID1, IID2, Z0, Z1, Z2 and PI_HAT absolutely required
  # POST: returns data.frame of all possible unique individuals' pairs, 
  #        whether they are related as specified in rel_pairs_detailed and
  #        the Z0,Z1,Z2 and PI_HAT
  
  possible_pairs <- permutePairs(individuals)
  rel <- logical(dim(possible_pairs)[1])
  values <- matrix(ncol=4,nrow=dim(possible_pairs)[1]) # to store corresponding Z0,Z1,Z2 and PI_HAT values
  phi_rel <- logical(dim(possible_pairs)[1])
  rel_type <- character(dim(possible_pairs)[1])
  
  for(i in 1:dim(possible_pairs)[1]) { # for each pair
    index<-which(as.character(rel_pairs_detailed[,"IID1"]) %in% possible_pairs[i,1])
    index2<-which(as.character(rel_pairs_detailed[,"IID2"]) %in% possible_pairs[i,1])
    if(length(index)>0) {
      for(j in 1:length(index)) {
        if(as.character(rel_pairs_detailed[index[j],"IID2"])==possible_pairs[i,2]) {
          rel[i]<-TRUE
          values[i,]<-as.numeric(rel_pairs_detailed[index[j],c("Z0","Z1","Z2","PI_HAT")])
          phi_rel[i] <- ifelse(as.numeric(values[i,4])>0.2,TRUE,FALSE)
          rel_type[i] <- getRelationship_type(as.numeric(values[i,1:3]),as.numeric(values[i,4]))
        }
      }
    }
    if(length(index2)>0) {
      for(k in 1:length(index2)) {
        if(as.character(rel_pairs_detailed[index2[k],"IID1"])==possible_pairs[i,2]) {
          rel[i]<-TRUE
          values[i,]<-as.numeric(rel_pairs_detailed[index2[k],c("Z0","Z1","Z2","PI_HAT")])
          phi_rel[i] <- ifelse(as.numeric(values[i,4])>0.2,TRUE,FALSE)
          rel_type[i] <- getRelationship_type(as.numeric(values[i,1:3]),as.numeric(values[i,4]))
        }
      }
    }    
  }
  result_table <- as.data.frame(cbind(possible_pairs,rel,values,phi_rel,rel_type))
  colnames(result_table)<-c("Pair1","Pair2","Pair_in_file","Z0","Z1","Z2", "PI_HAT","Related","Relationship_type")
  return(result_table)
}
#----------------------------------------------------------------------------------------------
relationship_summary2 <- function(individuals, rel_pairs_detailed){ #same but includes FIDs
  # PRE: individuals contains IIDs to be tested for relatedness as specified in rel_pairs_detailed;
  #      rel_pairs is the genome file PLINK output with IID1, IID2, Z0, Z1, Z2 and PI_HAT absolutely required
  # POST: returns data.frame of all possible unique individuals' pairs, 
  #        whether they are related as specified in rel_pairs_detailed and
  #        the Z0,Z1,Z2 and PI_HAT
  
  possible_pairs <- permutePairs(individuals)
  fids1 <- character(dim(possible_pairs)[1])
  fids2 <- character(dim(possible_pairs)[1])
  rel <- logical(dim(possible_pairs)[1])
  values <- matrix(ncol=4,nrow=dim(possible_pairs)[1]) # to store corresponding Z0,Z1,Z2 and PI_HAT values
  phi_rel <- logical(dim(possible_pairs)[1])
  rel_type <- character(dim(possible_pairs)[1])
  
  for(i in 1:dim(possible_pairs)[1]) { # for each pair
    index<-which(as.character(rel_pairs_detailed[,"IID1"]) %in% possible_pairs[i,1])
    index2<-which(as.character(rel_pairs_detailed[,"IID2"]) %in% possible_pairs[i,1])
    if(length(index)>0) {
      for(j in 1:length(index)) {
        if(as.character(rel_pairs_detailed[index[j],"IID2"])==possible_pairs[i,2]) {
          rel[i]<-TRUE
          fids1[i] <- as.character(rel_pairs_detailed[index[j],"FID1"])
          fids2[i] <- as.character(rel_pairs_detailed[index[j],"FID2"])
          values[i,]<-as.numeric(rel_pairs_detailed[index[j],c("Z0","Z1","Z2","PI_HAT")])
          phi_rel[i] <- ifelse(as.numeric(values[i,4])>0.2,TRUE,FALSE)
          rel_type[i] <- getRelationship_type(as.numeric(values[i,1:3]),as.numeric(values[i,4]))
        }
      }
    }
    if(length(index2)>0) {
      for(k in 1:length(index2)) {
        if(as.character(rel_pairs_detailed[index2[k],"IID1"])==possible_pairs[i,2]) {
          rel[i]<-TRUE
          fids1[i] <- as.character(rel_pairs_detailed[index2[k],"FID1"])
          fids2[i] <- as.character(rel_pairs_detailed[index2[k],"FID2"])
          values[i,]<-as.numeric(rel_pairs_detailed[index2[k],c("Z0","Z1","Z2","PI_HAT")])
          phi_rel[i] <- ifelse(as.numeric(values[i,4])>0.2,TRUE,FALSE)
          rel_type[i] <- getRelationship_type(as.numeric(values[i,1:3]),as.numeric(values[i,4]))
        }
      }
    }    
  }
  result_table <- as.data.frame(cbind(fids1,possible_pairs[,1],fids2,possible_pairs[,2],rel,values,phi_rel,rel_type))
  colnames(result_table)<-c("FID1","IID1", "FID2", "IID2","Pair_in_file","Z0","Z1","Z2", "PI_HAT","Related","Relationship_type")
  return(result_table)
}
#----------------------------------------------------------------------------------------------
removeRelated <- function(individuals, rel_pairs) {
  # PRE: individuals is a character vector of IIDs of length>1; 
  #      rel_pairs is a data.frame containing at least the following columns:
  #      IID1, IID2, Z0, Z1, Z2 and PI_HAT
  # POST: removes individuals that are related to each other as specified in rel_pairs
  #       and returns all non-related individuals. Returns NULL if all are related.
  result<-NULL
  if(allRelated(individuals,rel_pairs)) {
    return(NULL) #All provided individuals are related
  } else {
    rel_table <- relationship_summary(individuals, rel_pairs)
    rel_list <- subset(rel_table, rel_table$Related %in% TRUE)[,c("Pair1","Pair2")]
    rel_list <- unique(c(as.character(rel_list$Pair1),as.character(rel_list$Pair2)))
    non_rel_list <- subset(rel_table, rel_table$Related %in% FALSE)[,c("Pair1","Pair2")]
    non_rel_list<-unique(c(as.character(non_rel_list$Pair1),as.character(non_rel_list$Pair2)))
    result <- subset(non_rel_list,!(non_rel_list %in% rel_list))
  }
  if(length(result)==0) result<-NULL
  return(result)
}
#----------------------------------------------------------------------------------------------
getGrandparents <- function(inidividuals) {
  # PRE: individuals is a character vector of length>1
  #      the following global variables are available: sib_pairs, second_degree_pairs
  # POST: returns any grandparents in individuals
  #       IF NECESSARY, USES THE FACT THAT A "1" IN THE 6TH DIGIT OF THE IID INDICATES A GRANDPARENT
  
  no_sibs <- removeRelated(i_pairs,sib_pairs) # remove siblings
  if(length(no_sibs)==1) {
    if(as.numeric(substring(no_sibs,6,6)) < as.numeric(substring(related_individuals[i],6,6))) {
      return(no_sibs)
    } else if(as.numeric(substring(no_sibs,6,6)) >= as.numeric(substring(related_individuals[i],6,6))) {
      flagged[i]<-TRUE
      exceptions_count <- exceptions_count+1
      exceptions[i]<-exceptions_count
      print("\n")
      print(paste("Exception", exceptions_count,related_individuals[i],"flag3: ambiguous second-degree pair:",
                  related_individuals[i],no_sibs,
                  "sex:", sex[i],sex[which(related_individuals %in% no_sibs)]))
      print("\n")
      return(NULL)
    } 
  } else {
    no_second_degree <- removeRelated(i_pairs,rbind(sib_pairs,second_degree_pairs)) #naively remove 2nd degree relationships
    GPs <- NULL
    if(!is.null(no_second_degree)) { # both grandparents are present and unrelated
      return(no_second_degree)
    } else { # the case when a grandchild and one or two grandparents are present
      for(i in 1:length(no_sibs)) {
        if(as.numeric(substring(no_sibs[i],6,6))==1) GPs<-c(GPs,no_sibs[i])
      }
      return(GPs)
    } 
  }
}
#----------------------------------------------------------------------------------------------
assignParents <- function(ped, IID, parents) {
  # PRE: ped is the entire pedigree file of 6 columns
  #      IID is a character type vector of length 1 and is the individual for which we want to assign parents
  #      parents are the IIDs of the parents and can be a character type vector of length 1 or 2
  # POST: returns updated ped file with correspoinding father and mother IIDs for individual IID
  
  index <- which(as.character(ped[,2]) %in% IID)
  if(length(parents)==1){
    if(getSex(parents,fam)==1) { #Male
      ped[index,3]<-parents
    } else { #Female
      ped[index,4]<-parents 
    }
  } else { # two parents
    if(getSex(parents[1],fam)==1) {
      ped[index,3]<-parents[1]
      ped[index,4]<-parents[2]
    } else {
      ped[index,3]<-parents[2]
      ped[index,4]<-parents[1]
    }
  }
  return(ped)
}
#----------------------------------------------------------------------------------------------
assignParentsToSibPairs <- function(pedigree, all_sib_pairs, non_preferred_twin_list = NULL) {
  # PRE: ped has gone through assignment of existing parents and FIDs have been corrected
  # POST: to unify complex pedigrees, sib-pairs with missing parents will be assigned fake parent id's
  
  fake_count = 0
  newpedigree <- data.frame(FID=as.character(pedigree[,1]), IID=as.character(pedigree[,2]), PID=as.character(pedigree[,3]),
                       MID=as.character(pedigree[,4]), Sex=as.numeric(pedigree[,5]), Phenotype=as.numeric(pedigree[,6]), stringsAsFactors = F)
  
  sib1 <- as.character(all_sib_pairs[,'IID1'])
  sib2 <- as.character(all_sib_pairs[,'IID2'])
  if(!is.null(non_preferred_twin_list)) indexes <- which(sib1 %in% as.character(non_preferred_twin_list) | sib2 %in% as.character(non_preferred_twin_list))
  all_sib_pairs <- all_sib_pairs[-indexes,]
  
  for(i in 1:nrow(all_sib_pairs)) {  
    index1 <- which(newpedigree[,2] %in% as.character(all_sib_pairs[i,'IID1']))
    index2 <- which(newpedigree[,2] %in% as.character(all_sib_pairs[i,'IID2']))
    father1 <- newpedigree[index1,3]
    mother1 <- newpedigree[index1,4]
    father2 <- newpedigree[index2,3]
    mother2 <- newpedigree[index2,4]
    
    zero_count <- sum(c(father1 %in% "0", father2 %in% "0", mother1 %in% "0", mother2 %in% "0"))
    
    if(father1 %in% "0" & mother1 %in% "0" & father2 %in% "0" & mother2 %in% "0") {
      fake_count = fake_count + 1
      fakeid <- paste0("fake_parent",fake_count)
      newpedigree <- rbind(newpedigree, c(FID=newpedigree[index1,1], IID=fakeid, MID="0", PID="0", Sex=1, Phenotype=0))
      newpedigree[index1,3] <- fakeid
      newpedigree[index2,3] <- fakeid
      
      fake_count = fake_count + 1
      fakeid <- paste0("fake_parent",fake_count)
      newpedigree <- rbind(newpedigree, c(FID=newpedigree[index1,1], IID=fakeid, MID="0", PID="0", Sex=2, Phenotype=0))
      newpedigree[index1,4] <- fakeid
      newpedigree[index2,4] <- fakeid
    }
    
    else if( (father1 %in% "0" & father2 %in% "0") | (mother1 %in% "0" & mother2 %in% "0") ) {
      if(father1 %in% "0" & father2 %in% "0") {
        fake_count = fake_count + 1
        fakeid <- paste0("fake_parent",fake_count)
        newpedigree <- rbind(newpedigree, c(FID=newpedigree[index1,1], IID=fakeid, MID="0", PID="0", Sex=1, Phenotype=0))
        newpedigree[index1,3] <- fakeid
        newpedigree[index2,3] <- fakeid
      } else if(mother1 %in% "0" & mother2 %in% "0") {
        fake_count = fake_count + 1
        fakeid <- paste0("fake_parent",fake_count)
        newpedigree <- rbind(newpedigree, c(FID=newpedigree[index1,1], IID=fakeid, MID="0", PID="0", Sex=2, Phenotype=0))
        newpedigree[index1,4] <- fakeid
        newpedigree[index2,4] <- fakeid
      }
    }
    
    else if( ((father1 != "0" & mother1 != "0") | (father2 != "0" & mother2 != "0")) &
               zero_count != 0) {
      if(father1 != "0" & mother1 != "0") {
        newpedigree[index2,3] <- father1
        newpedigree[index2,4] <- mother1
      } else if(father2 != "0" & mother2 != "0") {
        newpedigree[index1,3] <- father2
        newpedigree[index1,4] <- mother2
      }
    }
    
    
    else if(zero_count != 0) {
      stop(paste("Problem with sibpair", all_sib_pairs[i,]))
    }
    
  }
  newpedigree2 <- data.frame(FID=newpedigree[,1], IID=newpedigree[,2], PID=newpedigree[,3], MID=newpedigree[,4],
                             Sex=as.numeric(newpedigree[,5]), Phenotype=as.numeric(newpedigree[,6]))
  return(newpedigree2)
}
#----------------------------------------------------------------------------------------------
assignFakeParentToDuos <- function(pedigree) {
  # PRE: a pedigree with missing parents assigned as "0"
  # POST: assigns a fake parent id to individuals who only have one parent (duos)
  
  fake_duo_count <- 0
  newpedigree <- data.frame(FID=as.character(pedigree[,1]), IID=as.character(pedigree[,2]), PID=as.character(pedigree[,3]),
                            MID=as.character(pedigree[,4]), Sex=as.numeric(pedigree[,5]), Phenotype=as.numeric(pedigree[,6]), stringsAsFactors = F)
    
  for(i in 1:nrow(newpedigree)) {
    father <- as.character(newpedigree[i,3])
    mother <- as.character(newpedigree[i,4])
    zero_count <- sum(c(father %in% "0", mother %in% "0"))
    if(zero_count == 1) { # we have a duo
      fake_duo_count = fake_duo_count + 1
      fakeid <- paste0("fakeduo", fake_duo_count)
      if(father %in% "0") {
        newpedigree <- rbind(newpedigree, c(FID=newpedigree[i,1], IID=fakeid, MID="0", PID="0", Sex=1, Phenotype=0))
        newpedigree[i,3] <- fakeid
      }
      else if(mother %in% "0") {
        newpedigree <- rbind(newpedigree, c(FID=newpedigree[i,1], IID=fakeid, MID="0", PID="0", Sex=2, Phenotype=0))
        newpedigree[i,4] <- fakeid
      }
    }
  }
  newpedigree2 <- data.frame(FID=newpedigree[,1], IID=newpedigree[,2], PID=newpedigree[,3], MID=newpedigree[,4],
                             Sex=as.numeric(newpedigree[,5]), Phenotype=as.numeric(newpedigree[,6]))
  return(newpedigree2)
}
############################ END OF FUNCTIONS #############################

################################# MAIN ####################################
related_individuals <- unique(c(as.character(parent_offsprings$IID1),as.character(parent_offsprings$IID2)))

#related_individuals <- subset(related_individuals, !(related_individuals %in% twins$IID1)) # remove one twin
#==============================================================================================================================
# THERE SHALL BE ONLY ONE SET OF PREFERRED TWIN FID/IID SET FOR THE 8 TWIN PAIRS
#==============================================================================================================================

preferred_twin <- c("1058_301", "7081_301", "7021_301", "1030_301", "8156_201", "8156_202", "4007_302", "1004_303")
preferred_twinFID <- c("1058", "7081", "7020", "1030", "8156", "8156", "4007", "1004")
non_preferred_twin <- c("1058_302", "7081_302", "7020_302", "1030_302", "7044_201", "7044_202", "4007_303", "1004_302")
non_preferred_twinFID <- c("1058", "7081", "7020", "1030", "7044", "7044", "4007", "1004")

related_individuals <- subset(related_individuals, !(related_individuals %in% non_preferred_twin)) # remove one twin

fids <- getFID(related_individuals, fam)
twins_index <- which(related_individuals %in% preferred_twin)


sex <- getSex(related_individuals, fam)
parent_offsprings <- subset(parent_offsprings, !(parent_offsprings$IID1 %in% non_preferred_twin | 
                                                   parent_offsprings$IID2 %in% non_preferred_twin)) # remove one twin
all_individuals <- unique(c(as.character(all_pairs$IID1),as.character(all_pairs$IID2)))
all_individuals <- subset(all_individuals, !(all_individuals %in% non_preferred_twin)) # remove one twin

theRest <- subset(all_individuals, !(all_individuals %in% related_individuals))
fidsRest <- getFID(theRest,fam)
sexRest <- getSex(theRest,fam)

is_parent <- logical(length(related_individuals))
flagged <- logical(length(related_individuals))
exceptions <- rep(0,length(related_individuals))
exceptions_count <- 0
exceptions_ipair <- rep(0,length(related_individuals))
exceptions_i2pair <- rep(0, length(related_individuals))
#sibs <- rbind(sib_pairs,second_degree_pairs)
sibs <- sib_pairs
sibs_and_second <- rbind(sib_pairs,second_degree_pairs)

ped <- data.frame(FID=c(fids,fidsRest),IID=c(related_individuals,theRest), PID=0,MID=0,Sex=c(sex,sexRest),phenotype=0) #all are parents by default
ped <- subset(ped, !(as.character(ped$IID) %in% as.character(non_preferred_twin)) ) # remove one twin

assigned_by_iid <- logical(dim(ped)[1])

# 1. 
# For all individuals that were identified to have parent-child relationships during Kinship Check (parent_offspring.txt file)
for (i in 1:length(related_individuals)) { # for each individual in parent-offspring relationship
  # find all related pairs for individual i
  i_pairs <- NULL
  i2_pairs <- NULL
  # 1.1
  if (!is_parent[i]) {
    # find all pairs related to i
    # 1.1.a
    i_pairs <- pairs(as.character(related_individuals[i]), 
                     data.frame(parent_offsprings$IID1, parent_offsprings$IID2))
    # 1.1.b
    # if only one individual paired with i
    if(length(i_pairs)==1) { # special case
      # 1.1.b.i.
      # find the relatives of the individual paired with i
      i2_pairs <- pairs(as.character(i_pairs), 
                        data.frame(parent_offsprings$IID1, parent_offsprings$IID2))
      # 1.1.b.ii
      if(length(i2_pairs)==1 & # if this individual also only pairs with one other individual
           related_individuals[i]==i2_pairs[1]) { # and is individual i
#         flagged[i]<-TRUE
#         flagged[which(related_individuals==i_pairs)]<-TRUE
#         exceptions_count <- exceptions_count+1
#         exceptions[i]<-exceptions_count
#         exceptions_ipair[which(related_individuals==i_pairs)]<-exceptions_count
#         print(paste("Exception", exceptions_count,"flag1: only one parent-child pair. related_individuals index:",i,
#                     "parent-child:",related_individuals[i],i_pairs,
#                     "sex:", sex[i],sex[which(related_individuals %in% i_pairs)]))
#         #stop()
        
        # 1.1.b.ii.1)
    # THIS CODE SPECIFIC TO RE PEDIGREES THANKS TO ID CODING SPECIFYING PARENT-CHILD RELATIONSHIP
        assigned_by_iid[i]<-TRUE
        i_code <- as.numeric(substring(related_individuals[i],6,6)) # get the 6th digit of IID eg. 3 for 7085_301
        pair_code <- as.numeric(substring(i_pairs,6,6))
        if(i_code < pair_code) {
          is_parent[i] <- TRUE
        } else if(i_code > pair_code) {
          ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
        
          # 1.1.b.ii.b.2)
        } else {
          flagged[i]<-TRUE
          #flagged[which(related_individuals==i_pairs)]<-TRUE
          exceptions_count <- exceptions_count+1
          exceptions[i]<-exceptions_count
          exceptions_ipair[which(related_individuals==i_pairs)]<-exceptions_count
          print("\n")
          print(paste("Exception", exceptions_count,related_individuals[i],"flag1: ambiguous parent-child pair:",
                      related_individuals[i],i_pairs,
                      "sex:", sex[i],sex[which(related_individuals %in% i_pairs)]))
          print("\n")
          #stop()
        }
        
      # 1.1.b.iii.
      } else if (length(i2_pairs)>1 &
                   allRelated(i2_pairs, data.frame(sibs$IID1,sibs$IID2))) { # case of when only one parent is enrolled
        is_parent[which(related_individuals==i_pairs)]<-TRUE # most likely that i_pairs is a parent
        ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
      
      # 1.1.b.iv. 
      } else if (length(i2_pairs)==2 &
                   !allRelated(i2_pairs, data.frame(sibs_and_second$IID1,sibs_and_second$IID2))) { # ie. i2_pairs are parents of i_pairs
        is_parent[i]<-TRUE
      
      # 1.1.b.v.
      } else { # We have grandfathers/grandmothers as well

#         flagged[i]<-TRUE
#         flagged[which(related_individuals==i_pairs)]<-TRUE
#         exceptions_count <- exceptions_count+1
#         exceptions[i]<-exceptions_count
#         exceptions_ipair[which(related_individuals==i_pairs)]<-exceptions_count
#         if(length(i2_pairs>1)) {
#           for (j in 1:length(i2_pairs)) {
#             flagged[which(related_individuals[i]==i2_pairs[j])]<-TRUE
#             exceptions_i2pair[which(related_individuals[i]==i2_pairs[j])]<-exceptions_count
#           }
#         }
        
        # 1.1.b.v.1)
    # THIS CODE SPECIFIC TO RE PEDIGREES THANKS TO ID CODING SPECIFYING PARENT-CHILD RELATIONSHIP
        assigned_by_iid[i]<-TRUE
        i_code <- as.numeric(substring(related_individuals[i],6,6)) # get the 6th digit of IID eg. 3 for 7085_301
        pair_code <- as.numeric(substring(i_pairs,6,6))
        if(i_code < pair_code) {
          is_parent[i] <- TRUE
        } else if(i_code > pair_code) {
          ped<-assignParents(ped,related_individuals[i],i_pairs)
        
        # 1.1.b.v.2)
        } else {
          flagged[i]<-TRUE
          #flagged[which(related_individuals==i_pairs)]<-TRUE
          exceptions_count <- exceptions_count+1
          exceptions[i]<-exceptions_count
          exceptions_ipair[which(related_individuals==i_pairs)]<-exceptions_count
          print("\n")
          print(paste("Exception", exceptions_count,related_individuals[i],"flag2: This case cannot be resolved.",sep=" "))
          print(paste("One parent-child pair:", related_individuals[i],i_pairs, "Sex:", sex[i],sex[which(related_individuals %in% i_pairs)],sep=" "))
          print(paste("and",i_pairs,"\' parent or child:",sep=" "))
          print(paste(i2_pairs)) 
          print(paste("Sex:",sex[which(related_individuals %in% i2_pairs)],sep=" "))
          print("Summary:")
          print(relationship_summary(i2_pairs,rbind(sibs,second_degree_pairs,unrelated2)))
          print("\n")
          #stop()
        }
        
        
      }
    
  # 2.a.
  # If more than one individual pairs with i
  } else if(length(i_pairs)>1) { # if more than one individual paired with i
    # 2.a.i.
      if(allRelated(i_pairs, data.frame(sibs$IID1,sibs$IID2))) { # if they are all sibs
        is_parent[i]<-TRUE # i is most likely a parent
      
      # 2.a.ii.
      } else if (length(i_pairs)==2 & # if there are two individuals related with i
                   !allRelated(i_pairs, data.frame(sibs_and_second$IID1,sibs_and_second$IID2)) & # and they are not related
                   sex[which(related_individuals==i_pairs[1])] != sex[which(related_individuals==i_pairs[2])]) { 
                    # and the paired individuals are of opposite sex
        ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
      
      # 2.a.iii.
      } else { # we may have grandparents
        # 2.a.iii.1)
        grandParents <- getGrandparents(i_pairs)
        
        if (is.null(grandParents) | length(grandParents)>2) { # most likely half-sibs or ambiguous pair
          is_parent[i]<-TRUE
        } else {
          ped<-assignParents(ped,related_individuals[i],as.character(grandParents))
        }
      }
    # 2.b.
    } else {
      stop(paste("No pairing for individual ",related_individuals[i]))
    }
  }
}

# SOLVE AND ASSIGN AMBIGUOUS PAIRS MANUALLY
# "Exception 1 7038_905 flag1: ambiguous parent-child pair: 7038_905 7038_901 sex: 1 2"
ped1 <- assignParents(ped, "7038_905", "7038_901")

# Exception 2 7055_202 flag3: ambiguous second-degree pair: 7055_202 7055_302 sex: 2 1
# Leave 7055_202 as is (parent of 7055_302 with 7055_201 as father)

# Exception 2 7033_904 flag1: ambiguous parent-child pair: 7033_904 7033_905 sex: 2 2
ped2 <- assignParents(ped1, "7033_904", "7033_905")

# Exception 3 7004_101 flag3: ambiguous second-degree pair: 7004_101 7004_203 sex: 2 2
# Leave 7004_101 as is (7004_101 is parent of 7004_203)

# Exception 3 7038_901 flag1: ambiguous parent-child pair: 7038_901 7038_905 sex: 2 1
# Leave as is as we have already assigned 7038_901 as the parent of 7038_905

# Exception 4 7033_905 flag1: ambiguous parent-child pair: 7033_905 7033_904 sex: 2 2"
# Leave as is as we have already assigned 7033_905 as the parent of 7033_904

newped <- ped2

#setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/")
#write.table(newped,"17_REped_IIDs_updated_twin_correction.txt", quote=FALSE,row.names=FALSE, col.names=TRUE)



#=====================================================================
# MUST UPDATE FIDs TO REFLECT RELATIONSHIPS IN THE SAME FAMILY
#=====================================================================
unmatched_parent_offsprings <- parent_offsprings[which(as.character(parent_offsprings$FID1) != as.character(parent_offsprings$FID2)),]
unmatched_sib_pairs <- sib_pairs[which(sib_pairs$FID1 != sib_pairs$FID2),]
unmatched_second_degree_pairs <- second_degree_pairs[which(second_degree_pairs$FID1 != second_degree_pairs$FID2),]

for(i in 1:nrow(unmatched_parent_offsprings)) {
  index <- which(as.character(newped[,2]) %in% unmatched_parent_offsprings[i,'IID2'])
  newped[index,1] <- as.character(unmatched_parent_offsprings[i,'FID1'])
}
for(i in 1:nrow(unmatched_sib_pairs)) {
  index <- which(as.character(newped[,2]) %in% unmatched_sib_pairs[i,'IID2'])
  newped[index,1] <- as.character(unmatched_sib_pairs[i,'FID1'])
}
for(i in 1:nrow(unmatched_second_degree_pairs)) {
  index <- which(as.character(newped[,2]) %in% unmatched_second_degree_pairs[i,'IID2'])
  newped[index,1] <- as.character(unmatched_second_degree_pairs[i,'FID1'])
}


##########################################################################################
# 16-APR-2015 
# - Discovered that FID 7012 and 7016 are really just one huge family
# - 7021_301 belongs to family 7021 NOT 7020
# - Should assign fake parents to duos as well
##########################################################################################
newped[which(as.character(newped[,1]) %in% "7016"),1] <- "7012"
newped[which(as.character(newped[,2]) %in% "7021_301"),1] <- "7021"

newped2 <- data.frame(FID=as.character(newped$FID), IID=as.character(newped$IID), 
                      PID=as.character(newped$PID), MID=as.character(newped$MID), 
                      Sex=newped$Sex, Phenotype=newped$phenotype, stringsAsFactors = F)

#=====================================================================
# MUST ASSIGN FAKE PARENT IDs TO UNIFY PEDIGREES
#=====================================================================
newped3 <- assignParentsToSibPairs(newped2, sib_pairs, non_preferred_twin)

# Next identify duos and assign a fake parent
newped3 <- assignFakeParentToDuos(newped3)



write.table(newped3, "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/REped_150416_update.txt", quote=F, row.names=F, col.names=F)


#===============================================================================================
library(kinship2)
#relation_matrix <- data.frame(id1=as.character(twins$IID1), id2=as.character(twins$IID2), code=rep(1,nrow(twins)))

buildRelationMatrix <- function(ped) {
  relation_matrix <- NULL
  spouse_pair <- NULL
  for(i in 1:nrow(ped)) {
    if(ped[i,'PID'] != "0" & ped[i,'MID'] != 0) {
      spouse_pair <- rbind(spouse_pair, 
                           data.frame(id1=ped[i,'PID'], id2=ped[i,'MID'], famid=ped[i,'FID'], code=4, 
                                      uniqid=paste0(ped[i,'PID'],"_",ped[i,'MID']), stringsAsFactors = F))
    }
  }
  relation_matrix <- spouse_pair[-which(duplicated(spouse_pair$uniqid)),c(1:4)]
  return(relation_matrix)
}


buildRelationMatrix2 <- function(ped) {
  relation_matrix <- NULL
  spouse_pair <- NULL
  ped[,'dadid'] <- ifelse(is.na(ped[,'dadid']),0,ped[,'dadid'])
  ped[,'momid'] <- ifelse(is.na(ped[,'momid']),0,ped[,'momid'])
  for(i in 1:nrow(ped)) {
    if(ped[i,'dadid'] != "0" & ped[i,'momid'] != 0) {
      spouse_pair <- rbind(spouse_pair, 
                           data.frame(id1=ped[i,'dadid'], id2=ped[i,'momid'], famid=ped[i,'fid'], code=4, 
                                      uniqid=paste0(ped[i,'dadid'],"_",ped[i,'momid']), stringsAsFactors = F))
    }
  }
  relation_matrix <- spouse_pair[-which(duplicated(spouse_pair$uniqid)),c(1:4)]
  return(relation_matrix)
}


# num_miss_fathers <- sum(is.na(newpid))
# num_miss_mothers <- sum(is.na(newmid))
# count_num_missing_parents <- sum(is.na(newpid)) + sum(is.na(newmid))
# fakeids <- NULL
# for(i in 1:count_num_missing_parents) fakeids <- c(fakeids, paste0("fakeid",i))
# count <- 1
# for(i in which(is.na(newpid))) {
#   newpid[i] <- fakeids[count]
#   count <- count + 1
# }
# for(i in which(is.na(newmid))) {
#   newmid[i] <- fakeids[count]
#   count <- count + 1
# }

# newfids <- c(as.character(newped3$FID), fakeids)
# newiids <- c(as.character(newped3$IID), fakeids)
# newsex <- c(newped3$Sex, rep(1,num_miss_fathers), rep(2, num_miss_mothers))
# newpid <- c(newpid, rep(NA, count_num_missing_parents))
# newmid <- c(newmid, rep(NA, count_num_missing_parents))

buildPedData <- function(ped) {
  newfids <- as.character(ped$FID)
  newiids <- as.character(ped$IID)
  newpid <- ifelse(as.character(ped$PID) == "0", NA, as.character(ped$PID))
  newmid <- ifelse(as.character(ped$MID) == "0", NA, as.character(ped$MID))
  newsex <- ped$Sex
  
  peddata <- data.frame(fid=newfids, id=newiids, 
                        dadid=newpid, momid=newmid, sex=newsex, status=rep(1,nrow(ped)))
  return(peddata)
}

buildPedData2 <- function(ped) {
  newfids <- as.character(ped$FID)
  newiids <- as.character(ped$IID)
  newpid <- ifelse(as.character(ped$PID) == "0", NA, as.character(ped$PID))
  newmid <- ifelse(as.character(ped$MID) == "0", NA, as.character(ped$MID))
  newsex <- ped$Sex
  
  peddata <- data.frame(fid=newfids, id=newiids, 
                        dadid=newpid, momid=newmid, sex=newsex, affected=ped$affected)
  
  return(peddata)
}

peddata <- buildPedData(newped3)
relation_matrix <- buildRelationMatrix(newped3)

ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=peddata$sex, 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
ped2 <- align.pedigree(ped)
#bitSize(ped)

# mfamid <- makefamid(id=peddata$fid, father.id=peddata$dadid, mother.id=peddata$momid)
# tmp <- familycheck(famid=newfids, id=newiids, father.id=newpid, mother.id=newmid, newfam=mfamid)
# subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
# subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree

# f <- attr(tmp, "join")['1001',]
# which(f != 0)

# newiids[mfamid==14]
# newiids[mfamid==15]
# newiids[mfamid==156]
# newiids[mfamid==157]



# Correcting family 1014 (2 individuals related at second-degree level)
newped4 <- data.frame(FID=as.character(newped3[,1]), IID=as.character(newped3[,2]), PID=as.character(newped3[,3]),
                      MID=as.character(newped3[,4]), Sex=as.numeric(newped3[,5]), Phenotype=as.numeric(newped3[,6]), stringsAsFactors = F)

newped4 <- rbind(newped4, c(FID="1014", IID="fakecustom1", PID="fakecustom3", MID="fakecustom4", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="1014", IID="fakecustom2", PID="0", MID="0", Sex=2, Phenotype=0))
newped4 <- rbind(newped4, c(FID="1014", IID="fakecustom3", PID="0", MID="0", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="1014", IID="fakecustom4", PID="0", MID="0", Sex=2, Phenotype=0))
newped4[which(newped4$IID %in% "1014_302"),3] <- "fakecustom1"
newped4[which(newped4$IID %in% "1014_302"),4] <- "fakecustom2"
newped4[which(newped4$IID %in% "1014_301"),3] <- "fakecustom3"
newped4[which(newped4$IID %in% "1014_301"),4] <- "fakecustom4"

relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree


# Correcting family 7023 (2 individuals who are actually unrelated)
newped4[which(newped4$IID %in% "7023_202"),1] <- "7023x"

relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree


# Correcting family 7033 (19 individuals with 4 splits)
# It is in fact 2 families
newped4 <- rbind(newped4, c(FID="7033", IID="fakecustom6", PID="0", MID="0", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7033", IID="fakecustom7", PID="0", MID="0", Sex=2, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7033", IID="fakecustom5", PID="fakecustom6", MID="fakecustom7", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7033", IID="fakecustom8", PID="0", MID="0", Sex=2, Phenotype=0))
newped4[which(newped4$IID %in% "7033_101"),3] <- "fakecustom6"
newped4[which(newped4$IID %in% "7033_101"),4] <- "fakecustom7"
newped4[which(newped4$IID %in% "7033_203"),3] <- "fakecustom5"
newped4[which(newped4$IID %in% "7033_203"),4] <- "fakecustom8"

newped4 <- rbind(newped4, c(FID="7033x", IID="fakecustom11", PID="0", MID="0", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7033x", IID="fakecustom12", PID="0", MID="0", Sex=2, Phenotype=0))
newped4[which(newped4$IID %in% "fake_parent27"),c(1,3)] <- c(FID="7033x", PID="fakecustom11")
newped4[which(newped4$IID %in% "fake_parent27"),c(1,4)] <- c(FID="7033x", MID="fakecustom12")
newped4[which(newped4$IID %in% "7033_905"),c(1,3)] <- c(FID="7033x", PID="fakecustom11")
newped4[which(newped4$IID %in% "7033_905"),c(1,4)] <- c(FID="7033x", MID="fakecustom12")
newped4[which(newped4$IID %in% c("7033_903", "fakeduo11", "7033_901", "7033_902", "7033_904")),1] <- "7033x"

relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree




# Correcting family 7038 (8 individuals with 2 splits and 1 unrelated)
# It is actually three families, with 7038_902 unrelated to the other two
newped4[which(newped4$IID %in% c("7038_901","fakeduo5","7038_905")),1] <- "7038x"
newped4[which(newped4$IID %in% c("7038_902")),1] <- "7038y"

relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree



# Correcting family 7085 (9 individuals with 2 splits)
newped4[which(newped4$IID %in% c("fakeduo43")),3] <- "7085_101"
newped4[which(newped4$IID %in% c("fakeduo43")),4] <- "fakeduo42"

relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree




subset(tmp, tmp$unrelated != 0 & tmp$split != 0)

# Correcting family 1001
# 1001_204 appears to be unrelated to the whole family
newped4[which(newped4$IID %in% c("1001_204")),1] <- "1001x"

# Correcting family 7039
newped4[which(newped4$IID %in% c("7039_903")),1] <- "7039x"
newped4[which(newped4$IID %in% c("7039_904")),1] <- "7039y"

# Correcting family 7040
newped4[which(newped4$IID %in% c("7040_905")),1] <- "7040x"

# Correcting family 7083
newped4[which(newped4$IID %in% c("7083_201")),1] <- "7083x"

# Correcting family 7086
newped4 <- rbind(newped4, c(FID="7086", IID="fakecustom13", PID="0", MID="0", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7086", IID="fakecustom14", PID="0", MID="0", Sex=2, Phenotype=0))
newped4[which(newped4$IID %in% c("fake_parent42")),3] <- "fakecustom13"
newped4[which(newped4$IID %in% c("fake_parent42")),4] <- "fakecustom14"
newped4[which(newped4$IID %in% c("7086_104")),3] <- "fakecustom13"
newped4[which(newped4$IID %in% c("7086_104")),4] <- "fakecustom14"


# Correcting family 7087
newped4 <- rbind(newped4, c(FID="7087x", IID="fakecustom15", PID="0", MID="0", Sex=1, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7087x", IID="fakecustom16", PID="0", MID="0", Sex=2, Phenotype=0))
newped4 <- rbind(newped4, c(FID="7087x", IID="fakecustom17", PID="0", MID="0", Sex=1, Phenotype=0))
newped4[which(newped4$IID %in% c("7087_202","7087_203")),1] <- "7087x"
newped4[which(newped4$IID %in% c("7087_202")),3] <- "fakecustom15"
newped4[which(newped4$IID %in% c("7087_202")),4] <- "fakecustom16"
newped4[which(newped4$IID %in% c("7087_203")),4] <- "fakecustom16"
newped4[which(newped4$IID %in% c("7087_203")),3] <- "fakecustom17"


relation_matrix <- buildRelationMatrix(newped4)
peddata <- buildPedData(newped4)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
tmp <- familycheck(famid=as.character(peddata$fid), id=peddata$id, father.id=peddata$dadid, mother.id=peddata$momid, newfam=mfamid)
subset(tmp, tmp$split !=1 & tmp$n != 1) # families with more than one split = no good; need to correct
subset(tmp, tmp$join != 0) # "join" will be non-zero if you are missing individuals in the family pedigree
subset(tmp, tmp$unrelated != 0 & tmp$split != 0)


hasDNA <- c(rep(TRUE,683),rep(FALSE,length(684:nrow(newped4))))
#pedigree.unrelated(ped, avail=hasDNA)


write.table(newped4, "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/REped_150416_update.txt", quote=F, row.names=F, col.names=F)



# Run RE_pheno_analysis2.R first
cts_set <- as.character(pheno_exome_cts$IID)
cts_ped_indexes <- which(newped4$IID %in% cts_set)

hasCTS <- rep(FALSE,nrow(newped4))
hasCTS[cts_ped_indexes] <- TRUE

#(cts_unrelated <- pedigree.unrelated(ped, avail=hasCTS))
# cts_unrelated2 <- cts_unrelated[1:117]
# temp <- relationship_summary2(cts_unrelated, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
# subset(temp, temp$Related %in% TRUE)

cts_or_ssd_indexes <- which(newped4$IID %in% pheno_exome_pass_cts_or_ssd$IID)
hasCTS_or_SSD <- rep(FALSE, nrow(newped4))
hasCTS_or_SSD[cts_or_ssd_indexes] <- TRUE
#(cts_or_ssd_unrelated <- pedigree.unrelated(ped, avail=hasCTS_or_SSD))
#cts_or_ssd_unrelated2 <- cts_or_ssd_unrelated[1:132]
#temp <- relationship_summary2(cts_or_ssd_unrelated2, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
#subset(temp, temp$Related %in% TRUE)


rdg_set <- as.character(pheno_exome_rdg$IID)
rdg_ped_indexes <- which(newped4$IID %in% rdg_set)
hasRDG <- rep(FALSE, nrow(newped4))
hasRDG[rdg_ped_indexes] <- TRUE
rdg_ped2 <- newped4[rdg_ped_indexes,]
#additional_ind <- 
# relation_matrix <- buildRelationMatrix(rdg_ped2)
# peddata <- buildPedData(rdg_ped2)
# ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
#                 status=peddata$status, relation=relation_matrix, famid = as.character(peddata$fid), missid="")
# mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
# 
#(rdg_unrelated <- pedigree.unrelated(ped, avail=hasRDG))
#rdg_unrelated2 <- rdg_unrelated[1:102]
#temp <- relationship_summary2(rdg_unrelated2, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
#subset(temp, temp$Related %in% TRUE)



cts_ped <- cbind(newped4[,1:5], status=c(rep(1,683),rep(0,length(684:nrow(newped4)))),
                 affected=as.integer(hasCTS))
cts_or_ssd_ped <- cbind(newped4[,1:5], status=c(rep(1,683),rep(0,length(684:nrow(newped4)))),
                        affected=as.integer(hasCTS_or_SSD))
rdg_ped <- cbind(newped4[,1:5], status=c(rep(1,683),rep(0,length(684:nrow(newped4)))),
                 affected=as.integer(hasRDG))

plotPedigree <- function(ped, famid) {
  # PRE: ped has 6 columns as in ped format
  # POST: famid is the family ID for which a plot is desired
  ped1 <- subset(ped, as.character(ped[,1]) %in% as.character(famid))
  peddata <- buildPedData2(ped1)
#  peddata$dadid <- ifelse(is.na(peddata$dadid),0,peddata$dadid)
#  peddata$momid <- ifelse(is.na(peddata$momid),0,peddata$momid)
#  newped1 <- data.frame(FID=ped1$FID, IID=peddata$id, PID=peddata$dadid, MID=peddata$momid, Sex=peddata$sex)
  relation_matrix <- buildRelationMatrix(ped1)
  ped2 <- pedigree(id=peddata$id, dadid=peddata$dadid, momid=peddata$momid, sex=as.numeric(peddata$sex), 
                  affected=peddata$affected)
  #ped2$id <- as.integer(gsub("\\d+\\_([0-9]*)","\\1",ped2$id))
  plot.pedigree(ped2,align=FALSE)
}

plotPedigree(rdg_ped, "1001")

# tmped <- buildPedData2(cts_or_ssd_ped)
# cts_or_ssd_ped2 <- subset(cts_or_ssd_ped, cts_or_ssd_ped$affected %in% 1)
# cts_or_ssd_pid <- unique(subset(cts_or_ssd_ped2$PID, !(cts_or_ssd_ped2$PID %in% c("0", cts_or_ssd_ped2$IID) )))
# cts_or_ssd_mid <- unique(subset(cts_or_ssd_ped2$MID, !(cts_or_ssd_ped2$MID %in% c("0", cts_or_ssd_ped2$IID) )))
# tmped <- buildPedData2(subset(cts_or_ssd_ped, cts_or_ssd_ped$affected %in% 1 | cts_or_ssd_ped$IID %in% c(cts_or_ssd_pid, cts_or_ssd_mid) ))
# kinmat <- kinship(id=tmped$id, dadid=tmped$dadid, momid=tmped$momid, sex=tmped$sex)

# for(i in 1:nrow(kinmat)) {
#   for(j in 1:ncol(kinmat)) {
#     if(i != j & kinmat[i,j] != 0) print(paste(rownames(kinmat)[i], colnames(kinmat)[j], kinmat[i,j]))
#   }
# }

#rdg_unrelated_ped <- subset(rdg_ped,rdg_ped$IID %in% rdg_unrelated2)
#cts_or_ssd_unrelated_ped <- subset(cts_or_ssd_ped, cts_or_ssd_ped$IID %in% cts_or_ssd_unrelated2)
#table(cts_or_ssd_unrelated_ped$affected)

relation_matrix <- buildRelationMatrix(rdg_ped)
peddata <- buildPedData(rdg_ped)
ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                affected=hasRDG)
mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))



unrel_list <- NULL
for(fam in unique(rdg_ped$FID)) { # for each of 201 families
  ped1 <- subset(rdg_ped, rdg_ped$FID %in% as.character(fam))
  
  if(sum(ped1$affected)==1) {
    unrel_list <- c(unrel_list, as.character(ped1[which(ped1$affected==1),2]))
  } 
  
  else if(nrow(ped1) != 1) {
    relation_matrix <- buildRelationMatrix(ped1)
    peddata <- buildPedData(ped1)
    ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                    affected=ped1$affected)
    mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
    
    genotypedVec <- ped1$status
    pheno <- ped1$affected
    #tmp <- findAvailAffected(ped, avail=genotypedVec, affstatus = ped1$affected)
    
    tmp <- pedigree.shrink(ped, avail=ped1$status, affected=ped1$affected)
    
    if(sum(tmp$pedObj$affected)!=0) {
      tmp2 <- pedigree.unrelated(tmp$pedObj, avail=tmp$pedObj$affected)
      unrel_list <- c(unrel_list, tmp2)
    } else if(length(tmp$pedObj$id)==0 & sum(ped1$affected)!=0) {
      print(fam)
    }
    
  } 
  else if(ped1$affected==1) {
    unrel_list <- c(unrel_list, as.character(ped1[,2]))
  }
}

# Manually add individuals in 7072
unrel_list <- c(unrel_list, "7072_201", "7072_202")
rdg_ped3 <- subset(rdg_ped, as.character(rdg_ped$IID) %in% as.character(unrel_list) )
temp <- relationship_summary2(unrel_list, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
subset(temp, temp$Related %in% TRUE)
table(rdg_ped3$affected)









unrel_list <- NULL
for(fam in unique(cts_or_ssd_ped$FID)) { # for each of 201 families
  ped1 <- subset(cts_or_ssd_ped, cts_or_ssd_ped$FID %in% as.character(fam))
  
  if(sum(ped1$affected)==1) {
    unrel_list <- c(unrel_list, as.character(ped1[which(ped1$affected==1),2]))
  } 
  
  else if(nrow(ped1) != 1) {
    relation_matrix <- buildRelationMatrix(ped1)
    peddata <- buildPedData(ped1)
    ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                    affected=ped1$affected)
    mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
    
    genotypedVec <- ped1$status
    pheno <- ped1$affected
    #tmp <- findAvailAffected(ped, avail=genotypedVec, affstatus = ped1$affected)
    
    tmp <- pedigree.shrink(ped, avail=ped1$status, affected=ped1$affected)
    
    if(sum(tmp$pedObj$affected)!=0) {
      tmp2 <- pedigree.unrelated(tmp$pedObj, avail=tmp$pedObj$affected)
      unrel_list <- c(unrel_list, tmp2)
    } else if(length(tmp$pedObj$id)==0 & sum(ped1$affected)!=0) {
      print(fam)
    }
    
  } 
  else if(ped1$affected==1) {
    unrel_list <- c(unrel_list, as.character(ped1[,2]))
  }
}

cts_or_ssd_ped3 <- subset(cts_or_ssd_ped, as.character(cts_or_ssd_ped$IID) %in% as.character(unrel_list) )
table(cts_or_ssd_ped3$affected)
temp <- relationship_summary2(unrel_list, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
subset(temp, temp$Related %in% TRUE)








unrel_list <- NULL
for(fam in unique(cts_ped$FID)) { # for each of 201 families
  ped1 <- subset(cts_ped, cts_ped$FID %in% as.character(fam))
  
  if(sum(ped1$affected)==1) {
    unrel_list <- c(unrel_list, as.character(ped1[which(ped1$affected==1),2]))
  } 
  
  else if(nrow(ped1) != 1) {
    relation_matrix <- buildRelationMatrix(ped1)
    peddata <- buildPedData(ped1)
    ped <- pedigree(id=as.character(peddata$id), dadid=as.character(peddata$dadid), momid=as.character(peddata$momid), sex=as.numeric(peddata$sex), 
                    affected=ped1$affected)
    mfamid <- makefamid(id=as.character(peddata$id), father.id=as.character(peddata$dadid), mother.id=as.character(peddata$momid))
    
    genotypedVec <- ped1$status
    pheno <- ped1$affected
    #tmp <- findAvailAffected(ped, avail=genotypedVec, affstatus = ped1$affected)
    
    tmp <- pedigree.shrink(ped, avail=ped1$status, affected=ped1$affected)
    
    if(sum(tmp$pedObj$affected)!=0) {
      tmp2 <- pedigree.unrelated(tmp$pedObj, avail=tmp$pedObj$affected)
      unrel_list <- c(unrel_list, tmp2)
    } else if(length(tmp$pedObj$id)==0 & sum(ped1$affected)!=0) {
      print(fam)
    }
    
  } 
  else if(ped1$affected==1) {
    unrel_list <- c(unrel_list, as.character(ped1[,2]))
  }
}

cts_ped3 <- subset(cts_ped, as.character(cts_ped$IID) %in% as.character(unrel_list) )
table(cts_ped3$affected)
temp <- relationship_summary2(unrel_list, rbind(sib_pairs, parent_offsprings, second_degree_pairs))
subset(temp, temp$Related %in% TRUE)


# write.table(newped4[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-FULL_pedigree_structure_specification.txt", quote=F, row.names=F, col.names=T)
# write.table(cts_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-CTS_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)
# write.table(cts_or_ssd_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-CTS_or_SSD_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)
# write.table(rdg_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-RDG_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)


ceu_individuals <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/CORRECTED_ANALYSIS/MANUAL_ETHNIC_ANALYSIS/35_ceu_tsi_set.txt",stringsAsFactors = F)[,2]  #561 individuals

ceu_rdg <- rdg_ped3[which(rdg_ped3[,2] %in% ceu_individuals),]
ceu_cts <- cts_ped3[which(cts_ped3[,2] %in% ceu_individuals),]
ceu_cts_ssd <- cts_or_ssd_ped3[which(cts_or_ssd_ped3[,2] %in% ceu_individuals),]

# write.table(cts_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-CTS_CEU_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)
# write.table(cts_or_ssd_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-CTS_or_SSD_CEU_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)
# write.table(rdg_ped3[,-6], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-RDG_CEU_best_unrelated_set.txt", quote=F, row.names=F, col.names=T)
# 



# My method
tmp <- relationship_summary2(as.character(rdg_ped2[,2]), rbind(sib_pairs,parent_offsprings, second_degree_pairs))
tmp <- subset(tmp, tmp$Related==TRUE)
rels_to_remove <- OptimizeRelatedRemoval(tmp, imiss)

rdg_unrel <- subset(rdg_ped2, !(rdg_ped2[,2] %in% rels_to_remove[,2]))
tmp <- relationship_summary2(as.character(rdg_unrel[,2]), rbind(sib_pairs,parent_offsprings, second_degree_pairs))
subset(tmp, tmp$Related==TRUE)

ceu_rdg_np_method <- subset(rdg_unrel, rdg_unrel[,2] %in% ceu_individuals)


#write.table(rdg_unrel, "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-RDG_CEU_unrelated_set_OptimizeRelatedInd_method.txt", quote=F, row.names=F, col.names=T)
