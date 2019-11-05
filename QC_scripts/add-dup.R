#!/usr/bin/R

args=(commandArgs(TRUE))

filename <- as.character(args[1])

rsupdate <- read.table(filename, header=F, stringsAsFactors=F)
rsupdate[which(duplicated(rsupdate[,1])),1] <- paste0(rsupdate[which(duplicated(rsupdate[,1])),1], ".dup")

write.table(rsupdate, filename, quote=F, row.names=F, col.names=F)


