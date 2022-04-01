#!/usr/bin/R
# This version simply runs flashPCA2 and saves the resulting eigenvalues and eigenvectors

args=(commandArgs(TRUE))

library(flashpcaR)
library(data.table)

bplink_filename <- as.character(args[1]) # root name without extension
out_prefix <- as.character(args[2]) # no extension needed

#fn <- file.path("data/intermediate_files/pca", bplink_filename)
print("Running flashPCA2")
f <- flashpca(bplink_filename, ndim=10)

print("Saving results")
ids <- rownames(f$vectors)
eigenvectors <- data.table(ids, f$vectors)
fwrite(as.data.frame(f$values), paste0(out_prefix,"_eigenvalues.txt"), 
       quote=F, row.names=T, col.names=F, sep="\t")
fwrite(eigenvectors, paste0(out_prefix,"_eigenvectors.txt"),
       quote=F, row.names=F, col.names=F, sep="\t")
fwrite(as.data.frame(f$pve), paste0(out_prefix,"_pve.txt"), 
       quote=F, row.names=T, col.names=F, sep="\t") # phenotypic variance explained
print("Done")
