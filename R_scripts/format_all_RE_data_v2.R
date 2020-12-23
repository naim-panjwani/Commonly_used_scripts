## 1. Extract all individuals genotyped on Exome and Omni

# Omni genotyped individuals:
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_original_data")
omni <- read.table("cnv_table_cut2.txt",header=FALSE)
omni_ind <- as.character(omni$V1) #285 individuals;

# Noticed some IIDs end in a "." and perhaps are duplicates...
index <- which(grepl("\\.",as.character(omni_ind)))
temp <- as.character(omni_ind[index])
temp2 <- gsub("(\\d*?)\\.","\\1", temp) #remove the dot
(length(which(as.character(temp2) %in% as.character(omni_ind))) == 0)
# TRUE, so no duplicates
# Replace with proper format:
omni_ind[index] <- temp2

# Exome chip genotyped individuals:
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v1")
exome_ind <- as.character(read.table("Step1_5ind_rm_UKRE_sampleID.fam",header=FALSE)$V2) #703 individuals; no duplicates


## 1.1 Mark each individual whether genotyped on Exome, Omni or both
genotyped_on_both <- subset(exome_ind, exome_ind %in% omni_ind)
exome_only_ind <- subset(exome_ind, !(exome_ind %in% genotyped_on_both))
omni_only_ind <- subset(omni_ind, !(omni_ind %in% genotyped_on_both))
all_genotyped_ind <- c(genotyped_on_both, exome_only_ind, omni_only_ind)

geno_table <- data.frame(IID=genotyped_on_both, Platform="Both")
geno_table <- rbind(geno_table, data.frame(IID=exome_only_ind, Platform="HumanExome"))
geno_table <- rbind(geno_table, data.frame(IID=omni_only_ind, Platform="HumanOmniExpress"))

## 1.2 Mark individuals that do not pass QC
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v4")
pass_QC_exome_ind <- as.character(read.table("Step7_twins_removed.fam", header=FALSE)$V2)
exome_ind_index <- which(geno_table$Platform %in% c("Both","HumanExome"))
exome_QC <- rep(NA, dim(geno_table)[1])
exome_QC[exome_ind_index] <- "FAIL"
exome_QC[as.character(geno_table$IID) %in% pass_QC_exome_ind] <- "PASS"

#######################  Need to do same after QC for Omni is done  ###############################
#setwd("")
#pass_QC_omni_ind <- as.character(read.table("",header=FALSE)$V2)
#omni_ind_index <- which(geno_table$Platform %in% c("Both","HumanOmniExpress"))
omni_QC <- rep(NA, dim(geno_table)[1])
#omni_QC[omni_ind_index] <- "FAIL"
#omni_QC[as.character(geno_table$IID) %in% pass_QC_omni_ind] <- "PASS"


geno_table <- cbind(geno_table, cbind(exome_QC, omni_QC))


## 1.3 Add ethnicity column to update later
ethnicity <- rep(NA, dim(geno_table)[1])


## 2. Extract and merge pheno data
setwd("/Users/naim panjwani/Documents/Strug/UKRE_original_data/")
pheno1 <- read.csv("GWAS-phenotypes10jul2014_modified.csv",header=TRUE)
pheno2 <- read.csv("Omni_without_phenotypes_mod.csv", header=TRUE)
pheno2 <- pheno2[,-2]
names(pheno1)[3] <- names(pheno2)[3]

pheno_common_ind <- subset(pheno1, as.character(pheno1$SAMPLE.ID) %in% as.character(pheno2$SAMPLE.ID))
# There is none
pheno_table <- rbind(pheno1, pheno2)
names(pheno_table)[1] <- "IID"


## 3. Identify individuals in pheno table whose genotype is missing and vice versa and mark them
missing_genotype <- subset(pheno_table, !(as.character(pheno_table$IID) %in% as.character(geno_table$IID)))
missing_genotype_ind <- as.character(missing_genotype$IID)
missing_genotype_index <- which(as.character(pheno_table$IID) %in% missing_genotype_ind)
genotyped <- rep(TRUE, dim(pheno_table)[1])
genotyped[missing_genotype_index] <- FALSE
pheno_table <- cbind(pheno_table, genotyped)

missing_phenotype <- subset(geno_table, !(as.character(geno_table$IID) %in% as.character(pheno_table$IID)))
missing_phenotype_ind <- as.character(missing_phenotype$IID)
missing_phenotype_index <- which(as.character(geno_table$IID) %in% missing_phenotype_ind)
phenotyped <- rep(TRUE, dim(geno_table)[1])
phenotyped[missing_phenotype_index] <- FALSE
geno_table <- cbind(geno_table, phenotyped)

# 3.1 Identify duplicated IIDs
needs.attention <- pheno_table[duplicated(as.character(pheno_table$IID)),]
needs.attention.index <- which(duplicated(as.character(pheno_table$IID)))
dup <- rep(FALSE, dim(pheno_table)[1])
dup[needs.attention.index] <- TRUE
pheno_table <- cbind(pheno_table, duplicated=dup)



## 4. Merge pheno and geno tables by IIDs
pheno_table <- pheno_table[,c(1,14,15,2:13)]
geno_table$IID <- as.character(geno_table$IID)
pheno_table$IID <- as.character(pheno_table$IID)
merged_table <- merge(geno_table, pheno_table, all=TRUE)
var_order <- c("IID", "Platform", "exome_QC","omni_QC","genotyped","phenotyped","duplicated","STUDY",
               "Ped.Singleton.Control","Site","prb","ageonset","cts","ssd","rdg","adhd","mig","aeds","aedpoly")
merged_table <- merged_table[,var_order]
merged_table$genotyped[is.na(merged_table$genotyped)] <- TRUE
merged_table$phenotyped[is.na(merged_table$phenotyped)] <- TRUE

write.csv(merged_table, "RE_merged_data.csv", quote=FALSE,row.names=FALSE)
