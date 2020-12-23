setwd("/Users/naim panjwani/Documents/Strug/R_code")

parent_offsprings <- read.table("parent_offspring.txt",header=TRUE)
sib_pairs <- read.table("siblings.txt", header=TRUE)
second_degree_pairs <- read.table("second_degree_relatives.txt",header=TRUE)
twins <- read.table("duplicates.txt",header=TRUE)
all_pairs <- read.table("Step14_kinship.genome",header=TRUE)

fam <- read.table("Step13_pruned_updatedFIDs.fam",header=FALSE)

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
  } else if(inRange(pi_hat,range(second_degree_pairs$PI_HAT))) {
    result<-"second_degree"
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
    return(no_sibs)
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

############################ END OF FUNCTIONS #############################

################################# MAIN ####################################
related_individuals <- unique(c(as.character(parent_offsprings$IID1),as.character(parent_offsprings$IID2)))
related_individuals <- subset(related_individuals, !(related_individuals %in% twins$IID2)) # remove one twin
fids <- getFID(related_individuals, fam)
sex <- getSex(related_individuals, fam)
parent_offsprings <- subset(parent_offsprings, !(parent_offsprings$IID1 %in% twins$IID2 | 
                                                   parent_offsprings$IID2 %in% twins$IID2)) # remove one twin
all_individuals <- unique(c(as.character(all_pairs$IID1),as.character(all_pairs$IID2)))
all_individuals <- subset(all_individuals, !(all_individuals %in% twins$IID2)) # remove one twin

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

ped <- data.frame(FID=c(fids,fidsRest),IID=c(related_individuals,theRest), PID=0,MID=0,Sex=c(sex,sexRest),phenotype=-999) #all are parents by default
ped <- subset(ped, !(as.character(ped$IID) %in% as.character(twins$IID2)) ) # remove one twin

for (i in 1:length(related_individuals)) { # for each individual in parent-offspring relationship
  # find all related pairs for individual i
  i_pairs <- NULL
  i2_pairs <- NULL
  if (!is_parent[i]) {
    # find all pairs related to i
    i_pairs <- pairs(as.character(related_individuals[i]), data.frame(parent_offsprings$IID1, parent_offsprings$IID2))
    
    # if only one individual paired with i
    if(length(i_pairs)==1) { # special case
      # find the relatives of the individual paired with i
      i2_pairs <- pairs(as.character(i_pairs), 
                        data.frame(parent_offsprings$IID1, parent_offsprings$IID2))
      if(length(i2_pairs)==1 & # if this individual also only pairs with one other individual (most likely i)
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
        
    # THIS CODE SPECIFIC TO RE PEDIGREES THANKS TO ID CODING SPECIFYING PARENT-CHILD RELATIONSHIP
        i_code <- as.numeric(substring(related_individuals[i],6,6)) # get the 6th digit of IID eg. 3 for 7085_301
        pair_code <- as.numeric(substring(i_pairs,6,6))
        if(i_code < pair_code) {
          is_parent[i] <- TRUE
        } else if(i_code > pair_code) {
          ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
        } else {
          flagged[i]<-TRUE
          flagged[which(related_individuals==i_pairs)]<-TRUE
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
        
      } else if (length(i2_pairs)>1 &
                   allRelated(i2_pairs, data.frame(sibs$IID1,sibs$IID2))) { # case of when only one parent is enrolled
        is_parent[which(related_individuals==i_pairs)]<-TRUE # most likely that i_pairs is a parent
        ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
      } else if (length(i2_pairs)==2 &
                   !allRelated(i2_pairs, data.frame(sibs_and_second$IID1,sibs_and_second$IID2))) { # ie. i2_pairs are parents of i_pairs
        is_parent[i]<-TRUE
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
        
    # THIS CODE SPECIFIC TO RE PEDIGREES THANKS TO ID CODING SPECIFYING PARENT-CHILD RELATIONSHIP
        i_code <- as.numeric(substring(related_individuals[i],6,6)) # get the 6th digit of IID eg. 3 for 7085_301
        pair_code <- as.numeric(substring(i_pairs,6,6))
        if(i_code < pair_code) {
          is_parent[i] <- TRUE
        } else if(i_code > pair_code) {
          ped<-assignParents(ped,related_individuals[i],i_pairs)
        } else {
          flagged[i]<-TRUE
          flagged[which(related_individuals==i_pairs)]<-TRUE
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
    } else if(length(i_pairs)>1) { # if more than one individual paired with i
      if(allRelated(i_pairs, data.frame(sibs$IID1,sibs$IID2))) { # if they are all sibs
        is_parent[i]<-TRUE # i is most likely a parent
      } else if (length(i_pairs)==2 & # if there are two individuals related with i
                   !allRelated(i_pairs, data.frame(sibs_and_second$IID1,sibs_and_second$IID2)) & # and they are not related
                   sex[which(related_individuals==i_pairs[1])] != sex[which(related_individuals==i_pairs[2])]) { 
                    # and the paired individuals are of opposite sex
        ped<-assignParents(ped,related_individuals[i],as.character(i_pairs))
        
      } else { # we may have grandparents
        grandParents <- getGrandparents(i_pairs)
        
        if (is.null(grandParents) | length(grandParents)>2) { # most likely half-sibs
          flagged[i]<-TRUE
          exceptions_count <- exceptions_count+1
          exceptions[i]<-exceptions_count
          for (j in 1:length(i_pairs)) {
            flagged[which(related_individuals==i_pairs[j])]<-TRUE
            exceptions_ipair[which(related_individuals==i_pairs[j])]<-exceptions_count
          }
          print("\n")
          print(paste("Exception", exceptions_count,related_individuals[i],"flag3: One of the following is a grandparent:",sep=" "))
          print(paste(removeRelated(i_pairs,sibs),sep=" "))
          print(paste("Parent-child relationship with:",related_individuals[i],sep=" "))
          print("Summary:")
          print(relationship_summary(c(related_individuals[i],i_pairs),rbind(sibs,second_degree_pairs,parent_offsprings,unrelated2)))
          print("\n")
          #stop()
        } else {
          ped<-assignParents(ped,related_individuals[i],as.character(grandParents))
        }
      }
    } else {
      stop(paste("No pairing for individual ",related_individuals[i]))
    }
  }
}

write.table(ped,"REped.txt", quote=FALSE,row.names=FALSE, col.names=TRUE)
