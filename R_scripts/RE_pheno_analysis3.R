# CAUTION
# Using tables() or xtabs()
# Some columns are character format in order to show NA's in the tables functions, but
# these NA's are character format (ie. "NA"), which is not the same as NA!!!



options(stringsAsFactors = F)
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data")

all_pheno_data <- read.table("RE_merged_data_correction2.csv",header=TRUE, sep=",", stringsAsFactors = F, na.strings = c("",NA))[-64,] # Get rid of one of the duplicated 7015_301 rows
all_pheno_data$cts <- ifelse((all_pheno_data$cts %in% c("","9","rt frontal")) | is.na(all_pheno_data$cts), "NA", all_pheno_data$cts)
all_pheno_data$ssd <- ifelse((all_pheno_data$ssd %in% c("","9")) | is.na(all_pheno_data$ssd), "NA", all_pheno_data$ssd)
all_pheno_data$rdg <- ifelse((all_pheno_data$rdg %in% c("","9")) | is.na(all_pheno_data$rdg), "NA", all_pheno_data$rdg)
all_pheno_data$adhd <- ifelse((all_pheno_data$adhd %in% c("","9")) | is.na(all_pheno_data$adhd), "NA", all_pheno_data$adhd)
all_pheno_data$mig <- ifelse((all_pheno_data$mig %in% c("","9")) | is.na(all_pheno_data$mig), "NA", all_pheno_data$mig)
all_pheno_data$aeds <- ifelse((all_pheno_data$aeds %in% c("","9")) | is.na(all_pheno_data$aeds), "NA", all_pheno_data$aeds)
all_pheno_data$aedpoly <- ifelse((all_pheno_data$aedpoly %in% c("","9")) | is.na(all_pheno_data$aedpoly), "NA", all_pheno_data$aedpoly)

all_pheno_data[grep("7015_201", all_pheno_data$IID), 'Site'] <- "USA"
# Since other family members in FID 7015 are from USA

#------------------------------------------------------------------------------------------------------------
applyIIDupdate <- function(fam, updatefile) {
  newfam <- fam
  famID <- paste(as.character(fam[,1]),as.character(fam[,2]),sep=":")
  oldFID <- paste(as.character(updatefile[,1]), as.character(updatefile[,2]),sep=":")
  index <- which(famID %in% oldFID)
  print(paste("Found", length(index), "of", length(oldFID), "changes to apply"))
  for(i in index) {
    updatefile_index <- which(oldFID %in% famID[i])
    newfam[i, 1] <- as.character(updatefile[updatefile_index,3])
    newfam[i, 2] <- as.character(updatefile[updatefile_index,4])
  }
  return(newfam)
}
#------------------------------------------------------------------------------------------------------------
getFID <- function(IID, fam) {
  fids <- character(length(IID))
  for(i in 1:length(IID)) {
    fids[i] <- as.character(fam[which(as.character(fam[,2]) %in%  as.character(IID)[i]),1])
  }
  return(fids)
}
#------------------------------------------------------------------------------------------------------------
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
  # Version 2
  # Initialize variables
  to_be_removed <- NULL
  fids <- NULL
  j <- 1
  
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
#------------------------------------------------------------------------------------------------
extract <- function(full_dataset, iids, compare_col_num) {
  resulting_table <- subset(full_dataset, as.character(full_dataset[,compare_col_num]) %in% as.character(iids))
  return(resulting_table)
}
#------------------------------------------------------------------------------------------------
setNAs <- function(x, nastrings = c("n/a","NA"), result="character") {
  x <- as.character(x)
  x <- ifelse(x %in% nastrings, NA, x)
  if(result == "numeric") {
    return(as.numeric(x))
  } else if(result == "character") {
    return(x)
  }
}
#------------------------------------------------------------------------------------------------



# HumanOmniExpress QC update 25-May-2015
omni_genotyped <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/UKRE_cnv_binary_TEST.fam", stringsAsFactors = F)
fam_omni_id_updates <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/14_id_updates.txt", stringsAsFactors = F)[-1,] # do not apply 3000x as want to merge with exome
fam_omni_id_updates2 <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/omni_id_updates2.txt", stringsAsFactors = F)
fam_omni_id_updates3 <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/omni_id_updates3.txt", stringsAsFactors = F)
omni_genotyped <- applyIIDupdate(omni_genotyped, fam_omni_id_updates)
omni_genotyped <- omni_genotyped[-which(omni_genotyped[,2] %in% "7042_302"),] # remove this duplicate as 7041_302 will be changed to 7042_302 for proper merging with exome
omni_genotyped <- applyIIDupdate(omni_genotyped, fam_omni_id_updates2)
omni_genotyped <- applyIIDupdate(omni_genotyped, fam_omni_id_updates3)
omni_pass_qc <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/16_UKRE_twins_removed.fam", stringsAsFactors = F)
omni_pass_qc <- applyIIDupdate(omni_pass_qc, fam_omni_id_updates)
omni_pass_qc <- applyIIDupdate(omni_pass_qc, fam_omni_id_updates2)
omni_pass_qc <- applyIIDupdate(omni_pass_qc, fam_omni_id_updates3)
omni_qc_fail <- subset(omni_genotyped, !(as.character(omni_genotyped[,2]) %in% as.character(omni_pass_qc[,2])))

# Are all Omni genotyped individuals in "all_pheno_data"?
subset(omni_genotyped, !(as.character(omni_genotyped[,2]) %in% as.character(all_pheno_data$IID)))
# RK025 is missing; it has high missingness on the Omni anyways though
# 7070x 7070_301x --> we don't know who this person could be as it is genotypically different in both chips (by pairwise IBS); believe to be mislabeled for Omni

all_pheno_data$omni_QC <- ifelse(as.character(all_pheno_data$IID) %in% as.character(omni_genotyped[,2]) & 
                                   as.character(all_pheno_data$IID) %in% as.character(omni_pass_qc[,2]), "PASS",
                                 ifelse(as.character(all_pheno_data$IID) %in% as.character(omni_genotyped[,2]) & 
                                          as.character(all_pheno_data$IID) %in% as.character(omni_qc_fail[,2]), "FAIL", NA))


# Take out individuals with rare epilepsies, lab controls and unknown phenotype individuals
pheno_data <- all_pheno_data[66:nrow(all_pheno_data),-7]
# ADD 7015_301 and 7015_201 back in
pheno_data <- rbind(all_pheno_data[62:63,-7],pheno_data)

fam_exome <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/Step1_5ind_rm_UKRE_sampleID.fam", stringsAsFactors = F)
fam_exome_id_updates <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/CORRECTED_ANALYSIS/MANUAL_ETHNIC_ANALYSIS/150420-CORRECTIONS/id_updates.txt", stringsAsFactors = F)
exome_genotyped <- applyIIDupdate(fam_exome, fam_exome_id_updates)


fam_merge <- merge(omni_genotyped[,1:2], exome_genotyped[,1:2], all=T)

final_pheno_data <- cbind(FID=getFID(IID = pheno_data$IID, fam = fam_merge), pheno_data)


# Add ethnicity info
omni_eth <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/omnichip/QC_v3/ethnicities.txt", header=T, stringsAsFactors = F)
omni_eth[grep("3000x",omni_eth$FID),1] <- "3000"
names(omni_eth)[3] <- "Omni_Ethnicity"
exome_eth <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v5/CORRECTED_ANALYSIS/MANUAL_ETHNIC_ANALYSIS/35_ceu_tsi_set.txt", stringsAsFactors = F)
exome_eth <- cbind(exome_eth, exome_eth="CEU")
names(exome_eth) <- c("FID","IID", "Exome_CEU")

omni_merge <- merge(final_pheno_data, omni_eth, all.x=T, sort=F)
exome_merge <- merge(omni_merge, exome_eth, all.x=T, sort=F)

final_pheno_data <- exome_merge
final_pheno_ceu <- subset(final_pheno_data, grepl("CEU", final_pheno_data$Omni_Ethnicity) | grepl("American", final_pheno_data$Omni_Ethnicity) |
                            grepl("CEU", final_pheno_data$Exome_CEU))

#write.table(final_pheno_data, "150605-final_pheno_file.txt", quote=F, row.names=F, col.names=T, sep="\t")
#write.table(final_pheno_ceu, "150605-final_pheno_file_ceu_subset.txt", quote=F, row.names=F, col.names=T, sep="\t")

#================================================================================================================
#                              Add RE Affectedness Status -- ignore; may not be correct -- obtained from Progeny9
#================================================================================================================

RE_affectedness <- read.csv("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150925-RE_affectedness_mod.csv", stringsAsFactors=F)
RE_affectedness <- RE_affectedness[-nrow(RE_affectedness),]
RE_affectedness[,1] <- gsub("CK ","CK",RE_affectedness[,1])
RE_affectedness[,1] <- gsub("RK ","RK",RE_affectedness[,1])
RE_affectedness[,2] <- ifelse(RE_affectedness[,2] == 9, NA, RE_affectedness[,2])


# Ck117 ID change to CK117
final_pheno_data[grep("Ck",final_pheno_data[,1]),1:2] <- data.frame(FID="CK117", IID="CK117")


final_pheno_data2 <- merge(final_pheno_data, RE_affectedness, by="IID", sort=F, all.x=T)
final_pheno_ceu <- subset(final_pheno_data2, grepl("CEU", final_pheno_data2$Omni_Ethnicity) | grepl("American", final_pheno_data2$Omni_Ethnicity) |
                            grepl("CEU", final_pheno_data2$Exome_CEU))


#write.table(final_pheno_data, "150928-final_pheno_file2.txt", quote=F, row.names=F, col.names=T, sep="\t")
#write.table(final_pheno_ceu, "150028-final_pheno_file_ceu_subset2.txt", quote=F, row.names=F, col.names=T, sep="\t")

final_pheno_data <- final_pheno_data2



#========================================================================================================================
#                              Subset CEU CTS individuals per Site
#========================================================================================================================
all_imputed_ceu <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/150529-RDG_chr1_Exome_Omni_EUR/merging/11_Exome_Omni_SS_rdg_dbSNP141_id_update3.fam", stringsAsFactors = F)
avail_ceu_RE <- subset(final_pheno_ceu, final_pheno_ceu$IID %in% all_imputed_ceu[,2])
ss_controls <- subset(all_imputed_ceu, all_imputed_ceu[,6] %in% 1)
names(ss_controls) <- c("FID","IID","MID","PID","Sex","cts")

setwd("/Users/naim/Documents/Strug/UKRE/UKRE_assoc/150528-CTS_ELP4_Exome_Omni_EUR/Site_analysis/")
# Canada+USA CEU CTS = 69 USA + 13 Canadians = 82
north_american_ceu_cts <- subset(avail_ceu_RE, (avail_ceu_RE$Site %in% c("USA", "Canada")) & 
                                   (avail_ceu_RE$cts %in% 1))
temp <- north_american_ceu_cts[,c('FID','IID','cts')]
temp[,'cts'] <- as.integer(temp[,'cts']) + 1
# write.table(rbind(temp, ss_controls[,c('FID','IID','cts')]), 
#             "150605-north_american_cts_ceu_plus_749SS.txt", quote=F, row.names=F, col.names=F)
 
# France+UK+Argentina CEU CTS = 7 French + 59 UK + 4 Argentina = 70
european_ceu_cts <- subset(avail_ceu_RE, (avail_ceu_RE$Site %in% c("UK", "France")) &
                             (avail_ceu_RE$cts %in% 1))
temp <- european_ceu_cts[,c('FID','IID','cts')]
temp[,'cts'] <- as.integer(temp[,'cts']) + 1
# write.table(rbind(temp, ss_controls[,c('FID','IID','cts')]), "150605-european_cts_ceu_plus_749SS.txt",
#             quote=F, row.names=F, col.names=F)

# Sardinians with CTS = 62
sardinians <- subset(avail_ceu_RE, (avail_ceu_RE$Site %in% "Sardinia") & (avail_ceu_RE$cts %in% 1))
temp <- sardinians[,c('FID','IID','cts')]
temp[,'cts'] <- as.integer(temp[,'cts']) + 1
# write.table(rbind(temp, ss_controls[,c('FID','IID','cts')]), "150605-sardinians_cts_plus_749SS.txt",
#             quote=F, row.names=F, col.names=F)



#================================================================================================================
#                              Subset CEU CTS|SSD individuals
#================================================================================================================

cts_or_ssd <- subset(avail_ceu_RE, avail_ceu_RE$cts %in% 1 | avail_ceu_RE$ssd %in% 1)
european_ceu_ctsssd <- subset(cts_or_ssd, cts_or_ssd$Site %in% c("UK", "France"))
north_american_ctsssd <- subset(cts_or_ssd, cts_or_ssd$Site %in% c("USA", "Canada"))
sardinian_ctsssd <- subset(cts_or_ssd, cts_or_ssd$Site %in% "Sardinia")

#================================================================================================================
#                              Non-proband with CTS (likely sib w/out RE)
#================================================================================================================
sib_w_no_RE <- subset(avail_ceu_RE, with(avail_ceu_RE, !(Site %in% "Sardinia") & cts %in% 1 & prb %in% 0))
write.table(sib_w_no_RE[,c(2,1,3:ncol(sib_w_no_RE))], "/Users/naim/Documents/Strug/UKRE/UKRE_original_data/150928-sib_pairs_without_RE.txt" ,
            row.names=F, col.names=T, quote=F)




#================================================================================================================
#                              Individuals in previous linkage studies
#================================================================================================================

elaine_data <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/old_2007_CTS_elp4_44_snps_fine_mapping/ElaineChr11SNP.txt", stringsAsFactors = F, header=F)
casecontrol_cts <- read.csv("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/old_2007_CTS_elp4_44_snps_fine_mapping/casecontrol_cts.csv", stringsAsFactors=F)
progeny_44snp_case_control_iids <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/old_2007_CTS_elp4_44_snps_fine_mapping/02_44_snp_finemapped_individuals.txt", stringsAsFactors=F, header=F)


# elaine_data[,2] <- gsub("C","CK",elaine_data[,2])  WRONG!!! CK-PREFIXED ID'S ARE KERALANS!!! ELAINE'S DATA ARE ALL CANADIANS
# elaine_data[,2] <- gsub("R","RK",elaine_data[,2]) WRONG!!! RK-PREFIXED ID'S ARE KERALANS!!!

cts_linkage_iids <- paste0(casecontrol_cts[,1], "_", casecontrol_cts[,2]) # ignore the 9000 individuals; they are independent controls
cts_ids <- progeny_44snp_case_control_iids[which(!(progeny_44snp_case_control_iids[,2] %in% cts_linkage_iids)),2]
cts_ids <- c(cts_ids, cts_linkage_iids)

ALL_IND_IN_PREV_STUDY <- c(elaine_data[,2], cts_ids)
cases_in_prev_study <- ALL_IND_IN_PREV_STUDY[grep("^R|^C|^\\d\\d\\d\\d_\\d\\d\\d$",ALL_IND_IN_PREV_STUDY)]

was_in_previous_study <- final_pheno_data[,1] %in% ALL_IND_IN_PREV_STUDY
final_pheno_data <- cbind(final_pheno_data, was_in_prev_study=was_in_previous_study)


# temp <- final_pheno_data[was_in_previous_study,]
# 148
# write.csv(xtabs(~Site+cts, data=temp), "~/Desktop/temp.csv", quote=F, row.names=T)


temp2 <- subset(cases_in_prev_study, !(cases_in_prev_study %in% final_pheno_data[,1]))
temp3 <- subset(temp2, temp2 %in% all_pheno_data[,1])
# really not anywhere to be found for these 175 individuals!

temp4 <- subset(temp2, !grepl("^C", temp2))


was_in_prev_study_ceu <- avail_ceu_RE[,1] %in% ALL_IND_IN_PREV_STUDY
avail_ceu_RE <- cbind(avail_ceu_RE, was_in_prev_study=was_in_prev_study_ceu)

# temp <- avail_ceu_RE[was_in_prev_study_ceu,]
# write.csv(xtabs(~Site+cts+prb, data=temp), "~/Desktop/temp.csv", quote=F, row.names=T)



#write.csv(xtabs(~Site+cts,data=final_pheno_data), "~/Desktop/temp.csv", quote=F, row.names=T)
#write.csv(xtabs(~Site+cts,data=avail_ceu_RE), "~/Desktop/temp.csv", quote=F, row.names=T)



# Overlap in association study
cts_or_ssd2 <- (subset(avail_ceu_RE, (avail_ceu_RE$cts %in% 1 | avail_ceu_RE$ssd %in% 1) & 
                         !(avail_ceu_RE$Site %in% "Sardinia")))
cts2 <- subset(avail_ceu_RE, avail_ceu_RE$cts %in% 1 & !(avail_ceu_RE$Site %in% "Sardinia"))
#186 cts or ssd cases non-sardinian
temp <- avail_ceu_RE[was_in_prev_study_ceu,]
length(temp2 <- intersect(temp[,1], cts_or_ssd2[,1]))
# 50 cts/ssd
temp3 <- subset(cts_or_ssd2, cts_or_ssd2[,1] %in% temp2)
length(unique(temp3[,2]))
#45 families cts/ssd
length(temp2 <- intersect(temp[,1], cts2[,1]))
# 47 cts
temp3 <- subset(cts2, cts2[,1] %in% temp2)
length(unique(temp3[,2]))
# 45 families cts


temp <- subset(cts_or_ssd2, !(cts_or_ssd2$Site %in% "Sardinia") &
                 cts_or_ssd2$Platform %in% c("HumanExome") )
table(temp$prb)
temp <- subset(cts_or_ssd2, !(cts_or_ssd2$Site %in% "Sardinia") &
                 cts_or_ssd2$Platform %in% c("HumanOmniExpress") )
table(temp$prb)
temp <- subset(cts_or_ssd2, !(cts_or_ssd2$Site %in% "Sardinia") &
                 cts_or_ssd2$Platform %in% c("Both") )
table(temp$prb)




# Lisa extracted the USA/UK site from 2007:
new_lisa_file <- read.csv("~/Documents/Strug/UKRE/UKRE_original_data/2007_study/phenotypes.csv", stringsAsFactors = F)
lisa_ids <- gsub("-","_", new_lisa_file[,3])
new_lisa_file[,3] <- lisa_ids

overlap <- subset(new_lisa_file, new_lisa_file[,3] %in% final_pheno_data[,1])
overlap_ceu <- subset(new_lisa_file, new_lisa_file[,3] %in% avail_ceu_RE[,1])


temp5 <- subset(cases_in_prev_study, !(cases_in_prev_study %in% overlap[,3]) & (cases_in_prev_study %in% final_pheno_data[,1]))
temp6 <- subset(final_pheno_data, final_pheno_data[,1] %in% temp5)
xtabs(~Site+cts+prb,data=temp6)


cts_2007 <- subset(new_lisa_file, new_lisa_file$CTS %in% 2)
temp <- cts_2007[which(cts_2007[,3] %in% final_pheno_data[,1]),]
temp2 <- cts_2007[which(!(cts_2007[,3] %in% temp[,3])),]
temp3 <- subset(final_pheno_data, final_pheno_data[,1] %in% temp2[,3])










# RDG
# liz_rdg <- read.table("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/Liz_RDG_linkage/pedin.7", stringsAsFactors = F, header=F)
# liz_rdg_iids <- liz_rdg[,2]
# indexes <- which(as.numeric(liz_rdg[,2]) < 10^6)
# liz_rdg_iids[indexes] <- as.numeric(liz_rdg_iids[indexes]*10)
# liz_rdg_iids <- gsub("([0-9][0-9][0-9][0-9])([0-9][0-9][0-9])","\\1_\\2", liz_rdg_iids)
# 
# was_in_RDG_study <- final_pheno_data[,1] %in% liz_rdg_iids
# #final_pheno_data[was_in_RDG_study,2]
# final_pheno_data <- cbind(final_pheno_data, was_in_RDG_study)
# 
# was_in_RDG_study_ceu <- avail_ceu_RE[,1] %in% liz_rdg_iids
# avail_ceu_RE <- cbind(avail_ceu_RE, was_in_RDG_study=was_in_RDG_study_ceu)





#================================================================================================================
#                              All rs662702 homozygous individuals in HumanCoreExome
#================================================================================================================
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data/")
rs662702_hom_list <- read.table("rs662702_homozygous_individuals.txt", header=T, stringsAsFactors = F)
rs662702_hom_list <- subset(rs662702_hom_list, rs662702_hom_list[,'rs662702_T'] %in% 2)
rs662702_cases_in_study <- read.csv("rs662702_homozygous_cases.txt", header=T, stringsAsFactors = F)

rs662702_hom_list <- applyIIDupdate(rs662702_hom_list, fam_exome_id_updates)

m <- merge(rs662702_hom_list, final_pheno_data, by="IID")
subset(rs662702_hom_list, !(rs662702_hom_list[,2] %in% m[,1]))

m2 <- merge(rs662702_hom_list, all_pheno_data, by="IID")

write.csv(m, "~/Desktop/rs662702_hom_RE.csv", quote=F, row.names=F)
write.csv(m2, "~/Desktop/rs662702_hom_rareRE.csv", quote=F, row.names=F)




#=================================================================================
#           Making a table for the paper
#=================================================================================

ethnicity <- ifelse(final_pheno_data$Omni_Ethnicity %in% "CEU" | final_pheno_data$Exome_CEU %in% "CEU",
                    "European", "Other")
final_pheno_data <- cbind(final_pheno_data, Ethnicity=ethnicity)

study_individuals <- subset(final_pheno_data, (final_pheno_data$Platform %in% c("Both", "HumanExome")) |
                              (final_pheno_data$Site %in% "Argentina" & final_pheno_data$cts %in% 1) |
                              final_pheno_data$IID %in% "7084_301")
#cts_assoc <- subset(study_individuals, study_individuals$cts %in% 1)
#xtabs(~Ethnicity, data=cts_assoc)
#cts_ssd_assoc <- subset(study_individuals, (study_individuals$cts %in% 1) | (study_individuals$ssd %in% 1))
#xtabs(~Ethnicity, data=cts_ssd_assoc)

# Merge sex
fam_merge2 <- merge(omni_genotyped, exome_genotyped, by="V2", all=T)
fam_merge2[847,'V5.x'] <- 2 # S226 Sardinian
fam_merge2$V5.x <- as.numeric(fam_merge2$V5.x)
sex <- fam_merge2$V5.x
for(i in 1:nrow(fam_merge2)) {
  if(is.na(fam_merge2$V5.x[i])) sex[i] <- fam_merge2$V5.y[i]
}
fam_merge2 <- cbind(fam_merge2, sex)
fam_merge2 <- fam_merge2[,c(1,12)]

names(fam_merge2) <- c("IID","Sex")
m <- merge(study_individuals, fam_merge2, by="IID",sort=F)

study_individuals <- m

temp <- ifelse(study_individuals$ageonset %in% "n/a", NA, study_individuals$ageonset)
study_individuals$ageonset <- temp

RE_computed <- NULL
hasRE <- function(prb, cts, ageonset, aeds, RE) {
  hasREresult <- NULL
  prb <- setNAs(prb, result="numeric")
  cts <- setNAs(cts, result="numeric")
  ageonset <- setNAs(ageonset, result="numeric")
  aeds <- setNAs(aeds, result="numeric")
  RE <- setNAs(prb, result="numeric")
  
  prb <- ifelse(is.na(prb), 0, prb)
  aeds <- ifelse(is.na(aeds), 0, aeds)
  RE <- ifelse(is.na(RE), 0, RE)
  if(prb == 1) {
    hasREresult <- TRUE
  } else if(prb == 0) {
    if(!is.na(ageonset)) {
      hasREresult <- TRUE
    } else if((cts == 1) & (aeds != 0)) {
      hasREresult <- TRUE
    } else if(RE == 1){
      hasREresult <- TRUE
    } else {
      hasREresult <- FALSE
    }
  }
  return(hasREresult)
}

for(i in 1:nrow(study_individuals)) {
  RE_computed <- c(RE_computed, hasRE(study_individuals$prb[i], 
                                      study_individuals$cts[i], 
                                      study_individuals$ageonset[i], 
                                      study_individuals$aeds[i], 
                                      study_individuals$RE[i]))
}

study_individuals <- cbind(study_individuals, RE_computed)

xtabs(~RE_computed+Sex, data=study_individuals)
xtabs(~RE_computed+cts, data=study_individuals)
xtabs(~RE_computed+ssd, data=study_individuals)
xtabs(~RE_computed+Ethnicity, data=study_individuals)
xtabs(~RE_computed+Site, data=study_individuals)





