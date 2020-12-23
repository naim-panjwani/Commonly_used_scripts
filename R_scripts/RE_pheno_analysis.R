setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data")

all_pheno_data <- read.csv("RE_merged_data_correction1.csv",header=TRUE)
pheno_data <- all_pheno_data[67:dim(all_pheno_data)[1],-7]
pheno_not_RE_study<-subset(pheno_data, !(grepl("RE.*",pheno_data$STUDY)))

# Phenotypes of those in the RE study
pheno_RE_study <- subset(pheno_data, grepl("RE.*",pheno_data$STUDY))

# Subset of those in HumanExome chip
pheno_exome <- subset(pheno_RE_study, pheno_RE_study$Platform %in% c("Both" ,"HumanExome"))

# Subset of those on HumanOmniExpress chip
pheno_omni <- subset(pheno_RE_study, pheno_RE_study$Platform %in% c("Both", "HumanOmniExpress"))

# Individuals per site
table(pheno_RE_study$Site)

# Do all probands have ageonset?
probands <- subset(pheno_RE_study, pheno_RE_study$prb %in% 1)


getFID <- function(IID, ind_data) {
  fids <- character(length(IID))
  for(i in 1:length(IID)) {
    fids[i] <- as.character(ind_data[which(as.character(IID)[i] == as.character(ind_data[,2])),1])
  }
  return(fids)
}

missing_summary <- function(x, category) {
  x <- as.character(x)
  x <- ifelse(x=="" | x=="9", NA, x)
  x <- as.factor(x)
  temp <- xtabs(~category+is.na(x))
  temp2 <- prop.table(table(category,is.na(x)),1)
  temp3 <- cbind(temp, temp2[,2]*100)
  colnames(temp3) <- c("Num complete","Num missing", "% Missing")
  totals <- c(sum(temp3[,1]), sum(temp3[,2]), (sum(temp3[,2])/(sum(temp3[,1])+sum(temp3[,2])))*100)
  temp3 <- rbind(temp3, Totals=totals)
  return(temp3)
}

detailed_missing_summary <- function(x, category) {  
  x <- as.character(x)
  x <- ifelse(x=="" | is.na(x), "NA", x)
  temp <- table(category, x)
#  names(attributes(temp)$dimnames) <- c("Site",varname)
  totals <- numeric(dim(temp)[2])
  percents <- numeric(dim(temp)[2])
  for (i in 1:dim(temp)[2]) {
    totals[i] <- sum(temp[,i])
  }
  for (i in 1:dim(temp)[2]) {
    percents[i] <- round((totals[i]/sum(totals))*100)
  }
  lastrow <- c(sum(percents[1:2]),NA,
               sum(percents[3:4]), NA)
  temp <- rbind(temp,Totals=totals, Percents=percents, Overall=lastrow)
  return(temp)
}


# Data completeness per site
# CTS
print(temp <- missing_summary(pheno_RE_study$cts, pheno_RE_study$Site), digits=3)
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_RE_study$cts, pheno_RE_study$Site))
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_exome$cts,pheno_exome$Site))
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_omni$cts,pheno_omni$Site))
write.csv(temp,"temp.csv",quote=F,row.names=T)


# SSD
print(temp <- missing_summary(pheno_RE_study$ssd, pheno_RE_study$Site), digits=3)
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_RE_study$ssd, pheno_RE_study$Site))
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_exome$ssd,pheno_exome$Site))
write.csv(temp,"temp.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_omni$ssd,pheno_omni$Site))
write.csv(temp,"temp2.csv",quote=F,row.names=T)


# RDG
print(temp <- missing_summary(pheno_RE_study$rdg, pheno_RE_study$Site), digits=3)
write.csv(temp,"temp2.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_RE_study$rdg, pheno_RE_study$Site))
write.csv(temp,"temp2.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_exome$rdg,pheno_exome$Site))
write.csv(temp,"temp2.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_omni$rdg,pheno_omni$Site))
write.csv(temp,"temp2.csv",quote=F,row.names=T)


# ADHD
print(temp <- missing_summary(pheno_RE_study$adhd, pheno_RE_study$Site), digits=3)
write.csv(temp,"temp2.csv",quote=F,row.names=T)

(temp<-detailed_missing_summary(pheno_RE_study$adhd, pheno_RE_study$Site))
(temp2<-detailed_missing_summary(pheno_exome$adhd,pheno_exome$Site))
(temp3<-detailed_missing_summary(pheno_omni$adhd,pheno_omni$Site))
(temp4<-rbind(temp,temp2,temp3))
write.csv(temp4,"temp4.csv",quote=F,row.names=T)

#=========================================================================
quick_summary <- function(pheno_col) {
  print(temp <- missing_summary(pheno_RE_study[,pheno_col], pheno_RE_study$Site), digits=3)
  write.csv(temp,"temp2.csv",quote=F,row.names=T)
  
  (temp<-detailed_missing_summary(pheno_RE_study[,pheno_col], pheno_RE_study$Site))
  (temp2<-detailed_missing_summary(pheno_exome[,pheno_col],pheno_exome$Site))
  (temp3<-detailed_missing_summary(pheno_omni[,pheno_col],pheno_omni$Site))
  (temp4<-rbind(temp,temp2,temp3))
  write.csv(temp4,"temp4.csv",quote=F,row.names=T)
}
#=========================================================================
quick_summary(which(names(pheno_RE_study) %in% "aedpoly"))
names(pheno_RE_study)


#============================================================================================================================
#===================================== Checking Assoc of Exome chip CTS individuals with ELP4 ===============================
#============================================================================================================================

# Those with CTS trait typed on the human exome and who have passed QC:
pheno_exome_cts <- subset(pheno_exome, pheno_exome$cts==1 & pheno_exome$exome_QC %in% "PASS") #174 individuals

# Isolate independent samples (unrelated individuals)
######## FIRST, Load ped_build_REversion2.R ########
#pheno_exome_cts_rels2 <- relationship_summary(as.character(pheno_exome_cts$IID), all_pairs) # takes too long to compare all pairs
pheno_exome_cts_rels <- relationship_summary2(as.character(pheno_exome_cts$IID), rbind(parent_offsprings,sibs_and_second))
pheno_exome_cts_rels <- subset(pheno_exome_cts_rels, pheno_exome_cts_rels$Related==TRUE)

#------------------------------------------------------------------------------------------------
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v3")
imiss <- read.table("Step2_missingness.imiss", header=TRUE)

#------------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------------
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
  to_be_removed <- NULL
  fids<-NULL
  j <- 1
  for (i in 1:dim(related_pairs)[1]) { # for each pair
    if(!(as.character(related_pairs$IID1[i]) %in% to_be_removed) & !(as.character(related_pairs$IID2[i]) %in% to_be_removed)) {
      current_FID <- as.character(related_pairs$FID1[i])
      if(current_FID != as.character(related_pairs$FID2[i])) stop("Relateds not in same family! FID1: ",related_pairs$FID1,"FID2: ",related_pairs$FID2)
      family_size <- countDuplicates(current_FID, as.character(related_pairs$FID1))
      if(countDuplicates(as.character(related_pairs$IID1[i]),c(as.character(related_pairs$IID1),as.character(related_pairs$IID2))) > 
                           countDuplicates(as.character(related_pairs$IID2[i]), c(as.character(related_pairs$IID1),as.character(related_pairs$IID2)))) {
        to_be_removed[j] <- as.character(related_pairs$IID1[i])
        fids[j]<-as.character(related_pairs$FID1[i])
        j<-j+1
      } else if (countDuplicates(as.character(related_pairs$IID1[i]),c(as.character(related_pairs$IID1),as.character(related_pairs$IID2))) == 
                   countDuplicates(as.character(related_pairs$IID2[i]), c(as.character(related_pairs$IID1),as.character(related_pairs$IID2)))) {
        temp <- WorstPair(related_pairs[i,],imiss)
        to_be_removed[j] <- temp[1,2]
        fids[j] <- temp[1,1]
        j<-j+1
      } else {
        to_be_removed[j] <- as.character(related_pairs$IID2[i])
        fids[j]<-as.character(related_pairs$FID2[i])
        j<-j+1
      }
    }
  }
  return(as.data.frame(cbind(FID=as.character(fids),IID=as.character(to_be_removed))))
}
#------------------------------------------------------------------------------------------------
rels_to_remove <- OptimizeRelatedRemoval(pheno_exome_cts_rels, imiss)
#write.csv(rels_to_remove, "temp.csv", quote=F, row.names=F)

# Remove the related individuals:
pheno_exome_cts_no_rels <- subset(pheno_exome_cts, !(pheno_exome_cts$IID %in% rels_to_remove[,2]))
# 151 individuals with CTS on the Exome chip and having no relatedness to each other

# Double-check that there are no related individuals:
relcheck <- relationship_summary2(as.character(pheno_exome_cts_no_rels$IID), rbind(parent_offsprings,sibs_and_second))
(relcheck <- subset(relcheck, relcheck$Related==TRUE))
# None found -- Good; it worked



setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v3")
fam <- read.table("Step1_5ind_rm_UKRE_sampleID.fam")

# IID corrections as of late:
fam[,2] <- as.character(fam[,2])
fam[which(as.character(fam[,2]) %in% "1035_2-1"),2] <- "1035_201"
fam[which(as.character(fam[,2]) %in% "1059_301dup"),2] <- "1059_302"
fam[which(as.character(fam[,2]) %in% "7015_301"),2] <- "7015_201"
fam[which(as.character(fam[,2]) %in% "7015_301dup"),2] <- "7015_301"
fam[which(as.character(fam[,2]) %in% "1061_302"),2] <- "1061_301"


fids<-getFID(pheno_exome_cts_no_rels$IID,fam)
pheno_exome_cts_no_rels <- cbind(FID=fids,pheno_exome_cts_no_rels)

# setwd("/Users/naim panjwani/Documents/Strug/UKRE_assoc//CTS_Exome_ELP4")
# write.table(pheno_exome_cts_no_rels[,1:2], "pheno_exome_cts_no_rels.txt",
#             quote=F, row.names=F, col.names=F)

# Check the ethnicity of these individuals:
setwd("/Users/naim panjwani/Documents/Strug/UKRE/UKRE_sampleID/pca_relatives_removed")
eth <- read.table("ethnicities.txt", header=T)

getEthnicity <- function(iid, eth) {
  iid <- as.character(iid)
  ethnicities <- NULL
  for (i in 1:length(iid)) {
    if(iid[i] %in% as.character(eth$IID)) {
      index <- which(as.character(eth$IID) %in% iid[i])
      ethnicities[i] <- as.character(eth$Ethnicity[index])
    } else {
      ethnicities[i] <- NA
    }
  }
  return(ethnicities)
}


head(pheno_exome_cts_no_rels2 <- cbind(pheno_exome_cts_no_rels, Ethnicity=getEthnicity(pheno_exome_cts_no_rels$IID, eth)))





#=============================================================================
extract <- function(full_dataset, iids, compare_col_num) {
  resulting_table <- subset(full_dataset, as.character(full_dataset[,compare_col_num]) %in% as.character(iids))
  return(resulting_table)
}
#=============================================================================


# Load identified CEU individuals among the Exome, CTS individuals (119 total)
setwd("/Users/naim panjwani/Documents/Strug/UKRE/UKRE_original_data")
ceu_cts_ind <- read.table("ceu_individuals.txt", header=F)

ceu_cts_ind_details <- extract(pheno_exome_cts_no_rels, ceu_cts_ind[,2], 2)
write.table(ceu_cts_ind_details, "ceu_cts_ind_details",
            quote=F, row.names=F, col.names=T)

temp<-table(ceu_cts_ind_details$Site)
temp<-ceu_cts_ind_details[ceu_cts_ind_details$prb==0,]
#write.csv(temp,"temp.csv",quote=F,row.names=F)

ceu_cts_UK <- subset(ceu_cts_ind_details, ceu_cts_ind_details$Site %in% "UK")
#write.table(ceu_cts_UK[,1:2], "ceu_cts_UK_ind.txt",
 #           quote=F, row.names=F, col.names=F)

ceu_cts_UKFrance <- subset(ceu_cts_ind_details, ceu_cts_ind_details$Site %in% c("UK","France") )
#write.table(ceu_cts_UKFrance[,1:2], "ceu_cts_UKFrance_ind.txt",
 #           quote=F, row.names=F, col.names=F)

ceu_cts_USA <- subset(ceu_cts_ind_details, ceu_cts_ind_details$Site %in% "USA")
#write.table(ceu_cts_USA[,1:2], "ceu_cts_USA_ind.txt",
 #           quote=F, row.names=F, col.names=F)

ceu_cts_excl_Canada <- subset(ceu_cts_ind_details, !(ceu_cts_ind_details$Site %in% "Canada"))
#write.table(ceu_cts_excl_Canada, "ceu_cts_excl_Canada_ind.txt",
 #           quote=F, row.names=F, col.names=F)

ceu_cts_France <- subset(ceu_cts_ind_details, ceu_cts_ind_details$Site %in% "France")
#write.table(ceu_cts_France, "ceu_cts_France_ind.txt",
 #          quote=F, row.names=F, col.names=F)

ceu_cts_Canada <- subset(ceu_cts_ind_details, ceu_cts_ind_details$Site %in% "Canada")
#write.table(ceu_cts_Canada, "ceu_cts_Canada_ind.txt",
 #          quote=F, row.names=F, col.names=F)

