#!/usr/bin/R
# This version simply runs PC-AiR and saves the resulting eigenvalues and eigenvectors

args=(commandArgs(TRUE))


bplink_filename <- as.character(args[1])
kinship_filename <- as.character(args[2])
out_prefix <- as.character(args[3])
# reference_panels_included <- as.logical(args[4])  # feature to be implemented in later versions

# source("https://bioconductor.org/biocLite.R")
# biocLite("GWASTools")
# biocLite("GENESIS")
# biocLite("SNPRelate")


library(GWASTools)
library(GENESIS)
library(SNPRelate)

#unlink("genotype.gds",force=T)


print("Converting binary PLINK file to GDS format")
snpgdsBED2GDS(bed.fn = paste0(bplink_filename,".bed"),
              bim.fn = paste0(bplink_filename,".bim"),
              fam.fn = paste0(bplink_filename,".fam"),
              out.gdsfn = paste0(bplink_filename,".gds"))
print("Reading genotype data")
geno <- GdsGenotypeReader(filename = paste0(bplink_filename,".gds"))
genoData <- GenotypeData(geno)

iids <- getScanID(genoData)
print("Reading KING kinship coefficients")
KINGmat <- king2mat(file.kin0 = paste0(kinship_filename,".kin0"),
                    file.kin = NULL, iids = iids)
# or run KING using snpgdsIBDKING() function

print("Running PC-AiR")
mypcair <- pcair(genoData = genoData, kinMat = KINGmat, divMat = KINGmat)

print("Generating PC-AiR plots for the first 3 PCs")
pdf(paste0(out_prefix,"_general_pca_plot_1of3.pdf"))
plot(mypcair, vx = 1, vy = 2)
dev.off()
pdf(paste0(out_prefix,"_general_pca_plot_2of3.pdf"))
plot(mypcair, vx = 2, vy = 3)
dev.off()
pdf(paste0(out_prefix,"_general_pca_plot_3of3.pdf"))
plot(mypcair, vx = 3, vy = 4)
dev.off()


eigenvalues <- mypcair$values
eigenvectors <- mypcair$vectors


str(mypcair)


print("Saving the eigenvalues and eigenvectors")
# Save eigenvalues and eigenvectors:
write.table(eigenvalues, paste0(out_prefix,"_eigenvalues.txt"), quote=F, row.names=T, col.names=F)
write.table(eigenvectors, paste0(out_prefix,"_eigenvectors.txt"), quote=F, row.names=T, col.names=F)
write.table(mypcair$rels, paste0(out_prefix,"_relateds.txt"), quote=F, row.names=F,col.names=F)
write.table(mypcair$unrels, paste0(out_prefix,"_unrelateds.txt"),quote=F,row.names=F,col.names=F)


print("Done")

# Later versions will enable the plotting of cases and controls (v2), and if included, the reference panels (v3)


