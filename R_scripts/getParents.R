#!/usr/bin/R

args=(commandArgs(TRUE))

famfile <- as.character(args[1])
full_famfile <- as.character(args[2])
use_integrated_full_fam <- as.logical(args[3])
dummy_parent_string <- as.character(args[4])
outfilename <- as.character(args[5])

if(is.na(use_integrated_full_fam)) use_integrated_full_fam <- FALSE
if(use_integrated_full_fam) {
  print("Using integrated full fam for UKRE study")
  #exome_fam <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150420-FULL_pedigree_structure_specification.txt", stringsAsFactors = F, header=T)
  #omni_fam <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/27_UKRE_id_updates3.fam", stringsAsFactors = F)
  exome_fam <- read.table("~/UKRE/UKRE_original_data/150420-FULL_pedigree_structure_specification.txt", stringsAsFactors = F, header=T)
  omni_fam <- read.table("~/UKRE/UKRE_QC/omnichip/QC_v3/27_UKRE_id_updates3.fam", stringsAsFactors = F)
  commons <- Reduce(intersect, list(exome_fam[,2], omni_fam[,2]))
  omni_minus_commons_fam <- subset(omni_fam, !(omni_fam[,2] %in% commons))[,1:5]
  names(omni_minus_commons_fam) <- names(exome_fam)
  full_fam <- rbind(exome_fam, omni_minus_commons_fam)  
} else {
  print(paste("Full fam file:", full_famfile))
  full_fam <- read.table(full_famfile, stringsAsFactors = F)
}

print(paste("Fam file:", famfile))
fam <- read.table(famfile, stringsAsFactors = F)


getFID <- function(IIDs, fam) {
  # PRE: IIDs is a vector of IID's contained within the fam file
  # POST: returns vector of FID's that matches the order of inputted IID's
  
  fids <- character(length(IIDs))
  for (i in 1:length(IIDs)) { # for each IID
    fids[i] <- as.character(fam[which(as.character(fam[,2])==as.character(IIDs[i])),1])
  }
  return(fids)
} # end of getFID


getParents <- function(fam, full_fam, dummy_parent_string = "fake") {
  # PRE: fam is a subset of full_fam columns 1-4 are FID, IID, Maternal_ID, Paternal_ID
  #      full_fam is the full set of individuals containing parental info
  # POST: gets the parents set of individuals specified in fam;
  #       will exclude individuals from the list of fam[,2], and those containing the dummy_parent_string
  
  parents <- unique(c(fam[,3], fam[,4]))
  parents <- parents[-which(parents %in% "0")]
  parents <- parents[-grep(dummy_parent_string, parents)]
  parents <- subset(parents, !(parents %in% fam[,2]))
  if(length(parents)==0) parents <- NULL
  return(parents)
}
# Example:
# parents <- getParents(subset(full_fam, full_fam[,2] %in% north_american_ceu_cts$IID), full_fam)


getAllParents <- function(fam, full_fam, dummy_parent_string = "fake") {
  # see getParents
  # This function recursively iterates to get parents AND 
  # parents of parents (ie. grandparents, great-grandparents, etc)
  
  parents <- getParents(fam, full_fam, dummy_parent_string)
  all_parents <- parents
  while(!is.null(parents)) {
    newfam <- subset(full_fam, full_fam[,2] %in% parents)
    parents <- getParents(newfam, full_fam, dummy_parent_string)
    parents <- subset(parents, !(parents %in% fam[,2] | parents %in% all_parents))
    all_parents <- c(all_parents, parents)
  }
  return(all_parents)
}
# Example:
# allparents <- getParents(subset(full_fam, full_fam[,2] %in% north_american_ceu_cts$IID), full_fam)


# Main
if(dummy_parent_string %in% "" | is.na(dummy_parent_string)) dummy_parent_string <- "fake"
print(paste("Dummy parents' string:", dummy_parent_string))
allparents <- getAllParents(fam, full_fam, dummy_parent_string)
allparentsfam <- cbind(getFID(allparents, full_fam), allparents)
if(is.null(allparents)) {
  print("No parents found")
} else {
  print(paste("Writing parents/grandparents list to", outfilename))
  write.table(allparentsfam, outfilename, quote=F, row.names=F, col.names=F)
}
