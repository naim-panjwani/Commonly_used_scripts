setwd("/Users/naim panjwani/Documents/Strug/UKRE_original_data/")

# Phenotype data
pheno <- read.csv("140624_phenotypes_copy.csv",header=TRUE)
pheno_rare <- subset(pheno, pheno$STUDY %in% c("Rare Epilepsies","MAE"))$SAMPLE.ID
pheno <- subset(pheno, pheno$STUDY %in% c("RE","Lab Control"))
progeny_pheno <- read.csv("RE_phenotypes.csv",header=TRUE)
progeny_pheno <- subset(progeny_pheno, 
                        !(as.character(progeny_pheno$Individual.name) %in% as.character(pheno_rare)))

# Cross-reference and check whether CTS values match between the two
progeny_common_ind <- subset(progeny_pheno, unique(as.character(progeny_pheno$Individual.name)) %in%
                               unique(as.character(pheno$SAMPLE.ID)))
pheno_common_ind <- subset(pheno, unique(as.character(pheno$SAMPLE.ID)) %in%
                             unique(as.character(progeny_pheno$Individual.name)))

progeny_common_ind <- cbind(as.character(progeny_common_ind$Individual.name), 
                            progeny_common_ind$CTS)
colnames(progeny_common_ind) <- c("IID","progeny_CTS")
progeny_common_ind <- progeny_common_ind[order(as.character(progeny_common_ind[,1])),]
pheno_common_ind <- cbind(as.character(pheno_common_ind$SAMPLE.ID),
                          pheno_common_ind$cts)
colnames(pheno_common_ind) <- c("IID", "xls_CTS")

progeny_common_ind <- as.data.frame(progeny_common_ind)
pheno_common_ind <- as.data.frame(pheno_common_ind)

comparison <- merge(progeny_common_ind,pheno_common_ind)
comparison[is.na(comparison$progeny_CTS),]
comparison[is.na(comparison$xls_CTS),]

differences <- comparison[ifelse(!is.na(comparison$progeny_CTS)!=!is.na(comparison$xls_CTS),1,0)==1,]
write.table(differences,"cts_differences_progeny.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(comparison,"cts_xls_progeny_comparison.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.csv(differences,"cts_differences_progeny.csv",
          quote=FALSE,row.names=FALSE)
write.csv(comparison,"cts_xls_progeny_comparison.csv",
          quote=FALSE,row.names=FALSE)


# Genotype data
cnv_study_table <- read.table("cnv_table_cut2.txt",header=FALSE)
cnv_ind <- cnv_study_table$V1
cnv_ind <- subset(cnv_ind, !(as.character(cnv_ind) %in% as.character(pheno_rare)))
setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3/")
exome_geno_individuals <- read.table("Step7_excludeMAF.fam",header=FALSE)$V2
exome_geno_individuals <- subset(exome_geno_individuals, !(as.character(exome_geno_individuals) %in%
                                                             as.character(pheno_rare)))

# Merge genotype data
common_ind <- subset(exome_geno_individuals, 
                     as.character(exome_geno_individuals) %in% as.character(cnv_ind))
cnv_minus_common <- subset(cnv_ind, 
                           !(as.character(cnv_ind) %in% as.character(common_ind)))
merged_ind <- c(as.character(exome_geno_individuals), as.character(cnv_minus_common))


# Merge phenotype data
pheno_comm_ind <- subset(as.character(progeny_pheno$Individual.name),
                         as.character(progeny_pheno$Individual.name) %in%
                           as.character(pheno$SAMPLE.ID))
#346 individuals in PROGENY already in phenotype sheet
pheno <- pheno[,c(1,9:15)]
progeny_pheno_reformat <- data.frame(SAMPLE.ID=progeny_pheno[,"Individual.name"],
                                     cts=progeny_pheno[,1],
                                     ssd=NA,
                                     rdg=progeny_pheno[,2],
                                     adhd=progeny_pheno[,3],
                                     mig=progeny_pheno[,5],
                                     aeds=NA,
                                     aedpoly=NA)
progeny_minus_common <- subset(progeny_pheno_reformat, 
                              !(as.character(progeny_pheno_reformat$SAMPLE.ID) %in% 
                                             as.character(pheno_comm_ind)))
sum(as.character(merged_ind) %in% as.character(progeny_minus_common))
# ZERO individuals for which the additional individuals from PROGENY coincide with available genotyped individuals
# ie. no need to merge with PROGENY data
merged_pheno <- rbind(pheno,progeny_minus_common)





# Compare pheno table with genotype data individuals
# 647 individuals in pheno file (excluding rare epilepsies/MAE)
# 871 genotyped individuals on file
# Exclude genotyped individuals with rare epilepsies/MAE
merged_ind <- subset(merged_ind, 
                           !(as.character(merged_ind) %in% as.character(pheno_rare)))
# Now have 816 individuals genotyped among those who do not have the rare epilepsies
# From these, how many are in the pheno file?
geno_and_pheno <- subset(merged_ind, 
                         as.character(merged_ind) %in% as.character(pheno$SAMPLE.ID))
# 632 out of the 816 genotyped individuals are present in the merged pheno file
missing_from_pheno_table <- subset(merged_ind, 
                                   !(as.character(merged_ind) %in% as.character(pheno$SAMPLE.ID)))
setwd("/Users/naim panjwani/Documents/Strug/UKRE_original_data/")
write.table(missing_from_pheno_table, "missing_complete_phenotypes.txt",
            quote=FALSE,row.names=FALSE, col.names=FALSE)

no_geno <- subset(pheno$SAMPLE.ID, 
                  !(as.character(pheno$SAMPLE.ID) %in% as.character(merged_ind)))

# Have genotype data on both chips for the following individuals:
cnv_ind[as.character(cnv_ind) %in% as.character(geno_and_pheno)]
sum(as.character(cnv_ind) %in% as.character(geno_and_pheno))
# 95 individuals



pheno_columns <- pheno[,c("cts","ssd","rdg","adhd","mig","aeds","aedpoly")]

all_pheno_missing <- as.character(pheno$SAMPLE.ID[which(apply(pheno_columns,1,function(x) sum(is.na(x) | x==9)==dim(pheno_columns)[2]))])
cts_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$cts) | pheno$cts==9])
ssd_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$ssd) | pheno$ssd==9])
rdg_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$rdg) | pheno$rdg==9])
adhd_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$adhd) | pheno$adhd==9])
mig_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$mig) | pheno$mig==9])
aeds_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$aeds) | pheno$aeds==9])
aedpoly_missing_list <- as.character(pheno$SAMPLE.ID[is.na(pheno$aedpoly) | pheno$aedpoly==9])

missing_numbers <- data.frame(num_individuals=length(merged_ind),all_pheno=length(all_pheno_missing), cts=length(cts_missing_list), ssd=length(ssd_missing_list), rdg=length(rdg_missing_list),
             adhd=length(adhd_missing_list), mig=length(mig_missing_list), aeds=length(aeds_missing_list), aedpoly=length(aedpoly_missing_list))
missing_ids <- list(all_pheno=all_pheno_missing, cts=cts_missing_list, ssd=ssd_missing_list, rdg=rdg_missing_list,
                    adhd=adhd_missing_list, mig=mig_missing_list, aeds=aeds_missing_list, aedpoly=aedpoly_missing_list)

tmp.wid<-getOption("width")
options(width=10000)
sink("missing_phenotypes_IIDs.txt")
print(missing_ids)
sink()
options(width=tmp.wid)



#### Which individuals on Omni (cnv study) have phenotype and which do not? 
cnv_study_table <- subset(cnv_study_table, !(as.character(cnv_study_table$V1) %in%
                                               as.character(pheno_rare)))
cnv_with_pheno_mergedset <- subset(cnv_study_table, as.character(cnv_study_table$V1) 
                                %in% as.character(merged_pheno$SAMPLE.ID) )
length(unique(as.character(merged_pheno$SAMPLE.ID)))
which(duplicated(as.character(merged_pheno$SAMPLE.ID)))
# 1 duplicate
which(duplicated(cnv_with_pheno_merged$V1))
# no duplicates
cnv_pheno <- subset(merged_pheno, as.character(merged_pheno$SAMPLE.ID) %in% 
                      cnv_with_pheno_mergedset$V1)
write.csv(cnv_pheno, "Omni_pheno.csv", 
          quote=FALSE, row.names=FALSE)

cnv_without_pheno_mergedset <- subset(cnv_study_table, !(as.character(cnv_study_table$V1) 
                                %in% as.character(merged_pheno$SAMPLE.ID)) )
write.csv(cnv_without_pheno_mergedset, "Omni_without_phenotypes.csv",
          quote=FALSE, row.names=FALSE)
