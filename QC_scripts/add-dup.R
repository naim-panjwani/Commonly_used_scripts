#!/usr/bin/R

library(data.table)

args=(commandArgs(TRUE))

filename <- as.character(args[1])  # bim file format

rsupdate <- fread(filename, header=F, stringsAsFactors=F)
rsupdate[which(duplicated(rsupdate[,2])),2] <- sapply(rsupdate[which(duplicated(rsupdate[,2])),2], function(x) paste0(x, ".dup"))

#newbim <- rsupdate[c(sapply(rsupdate[,2], function(x) !grepl(".dup",x)))]
#altsnpname <- sapply(1, function(i) paste0(newbim[,V1], ":", newbim[,V4], ":", newbim[,V5], ":", newbim[,V6]))[,1]

#setkey(newbim, V2)
#dups_i <- c(which(duplicated(altsnpname)), which(duplicated(altsnpname, fromLast=T)))
#dups2 <- NULL
#if(length(dups_i)>0) dups2 <- newbim[dups_i][,V2]


fwrite(rsupdate, filename, quote=F, row.names=F, col.names=F, sep="\t")


