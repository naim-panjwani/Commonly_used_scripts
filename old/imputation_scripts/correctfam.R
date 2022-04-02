#!/usr/bin/R

args=(commandArgs(TRUE))

fam_name <- as.character(args[1])
badfam_name <- as.character(args[2])

fam <- read.table(fam_name, stringsAsFactors=F)
badfam <- read.table(badfam_name, stringsAsFactors=F)

correct_fids <- NULL
correct_sex <- NULL
correct_mid <- NULL
correct_pid <- NULL

for(i in 1:nrow(badfam)) {
  famrow <- fam[which(fam[,2] %in% badfam[i,2]), ]
  correct_fids <- rbind(correct_fids, famrow[1])
  correct_sex <- rbind(correct_sex, famrow[5])
  correct_mid <- rbind(correct_mid, famrow[3])
  correct_pid <- rbind(correct_pid, famrow[4])
}

id_update <- cbind(oldFID=badfam[,1], oldIID=badfam[,2], newFID=correct_fids,newIID= badfam[,2])
sex_update <- cbind(FID=correct_fids, IID=badfam[,2], SEX=correct_sex)
parents_update <- cbind(FID=correct_fids, IID=badfam[,2], MID=correct_mid, PID=correct_pid)

write.table(id_update, "update_id.txt", quote=F, row.names=F, col.names=F)
write.table(sex_update, "update_sex.txt", quote=F, row.names=F, col.names=F)
write.table(parents_update, "update_parents.txt", quote=F, row.names=F, col.names=F)

