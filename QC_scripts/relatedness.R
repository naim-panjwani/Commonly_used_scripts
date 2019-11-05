#!/usr/bin/R

args=(commandArgs(TRUE))

relatednessfile <- as.character(args[1])
prefix <- as.character(args[2])
if(is.na(prefix)) {
  prefix <- ""
} else {
  prefix <- paste0(prefix,"_")
}

relatedness <- read.table(relatednessfile, header=TRUE, stringsAsFactors=F)
unrelated2 <- subset(relatedness, with(relatedness, Z2<0.23 & PI_HAT < 0.24))
unrelated2_pairs <- paste(as.character(unrelated2$IID1), as.character(unrelated2$IID2),sep=",")
second_degree_relatives <- subset(relatedness, with(relatedness, Z0>0.3 & Z2<0.1 & 
                                                      PI_HAT >=0.24 & PI_HAT <0.4))
(second_degree_pairs <- paste(as.character(second_degree_relatives$IID1), as.character(second_degree_relatives$IID2),sep=","))
sib_pairs <- subset(relatedness, with(relatedness, Z0>0.1 & Z1>0.3 & Z2>=0.1))
(sib_pairs_pairs <- paste(as.character(sib_pairs$IID1), as.character(sib_pairs$IID2),sep=","))
parent_offspring <- subset(relatedness, with(relatedness, Z0<0.1 & Z1>0.8 & Z2<0.2))
(parent_offspring_pairs <- paste(as.character(parent_offspring$IID1), as.character(parent_offspring$IID2),sep=","))
twins_duplicates <- subset(relatedness, with(relatedness, Z2>0.9))
(duplicate_pairs <- paste(as.character(twins_duplicates$IID1), as.character(twins_duplicates$IID2),sep=","))
(all_included <- dim(relatedness)[1] == (dim(unrelated2)[1]+dim(second_degree_relatives)[1]+dim(sib_pairs)[1]+
                                           dim(parent_offspring)[1]+dim(twins_duplicates)[1]) )

write.csv(second_degree_relatives, paste0(prefix, "second_degree_relatives.csv"),
            quote=FALSE,row.names=FALSE)
write.csv(sib_pairs, paste0(prefix, "siblings.csv"),
            quote=FALSE,row.names=FALSE)
write.csv(parent_offspring, paste0(prefix, "parent_offspring.csv"),
            quote=FALSE,row.names=FALSE)
write.csv(twins_duplicates, paste0(prefix, "duplicates.csv"),
            quote=FALSE,row.names=FALSE)
write.csv(unrelated2, paste0(prefix, "unrelated.csv"),
            quote=FALSE,row.names=FALSE)





