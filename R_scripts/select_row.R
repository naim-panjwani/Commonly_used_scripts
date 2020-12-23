# R script
args=(commandArgs(TRUE))

filename <- as.character(args[1])

data<-read.table(filename,header=F)
position <- data[nrow(data),1]
bpposition <- as.numeric(strsplit(as.character(position),"-")[[1]][2])
thisrow <- logical(nrow(data))
thisrow[length(thisrow)] <- TRUE
for(i in 1:(nrow(data)-1)) {
  if(((as.numeric(data[i,3]) - as.numeric(data[i,2])) == 1) & (as.numeric(data[i,3]) == as.numeric(bpposition))) {
    thisrow[i] <- TRUE
  }
}
result <- data[which(thisrow)[1],]
write.table(result, filename, quote=F, row.names=F, col.names=F)
