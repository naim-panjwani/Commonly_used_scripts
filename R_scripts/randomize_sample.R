setwd("/Users/naim panjwani/Desktop/UNC_CF_Ctrls")

fam <- read.table("06_UNC_CF_CTS_cases_CFLD_rm.fam",header=F)
fam_cases<-subset(fam,fam$V6 %in% 2)
fam <- subset(fam, fam$V6 %in% 1)
fam_males <- subset(fam, fam$V5 %in% 1)
fam_females<-subset(fam, fam$V5 %in% 2)
ids <- paste(fam$V1,fam$V2,sep=":")
ids_males <- paste(fam_males$V1,fam_males$V2,sep=":")
ids_females <- paste(fam_females$V1,fam_females$V2,sep=":")
ids_cases <- paste(fam_cases$V1, fam_cases$V2,sep=":")

randomized_males <- sample(ids_males)[1:363]
randomized_females <- sample(ids_females)[1:237]
randomized_sample <- c(randomized_males, randomized_females, ids_cases)
fids <- gsub("(.*):.*","\\1",randomized_sample)
iids <- gsub(".*:(.*)","\\1",randomized_sample)

final_sample <- cbind(fids,iids)
write.table(final_sample, "randomized_ctrls600_AND_all_cases.txt",
            quote=F, row.names=F, col.names=F)
