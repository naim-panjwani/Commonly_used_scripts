# Rscript to calculate best number of splits for pruning given the number of desired number of cores to use

args=(commandArgs(TRUE))

i <- c(1:100)
lhs <- (i^2+i)/2
rhs <- as.numeric(args[1])
num_splits<-length(which(lhs<rhs))
write.table(num_splits,"num_splits.txt",quote=F,row.names=F,col.names=F)
