# Inputs: observed association; permuted association; promoter or enhancer indicator matrix
# Outputs: contrast test statistic with empirical p-value

#args=(commandArgs(TRUE))

# permutation_filename <- as.character(args[1])
# observed_assoc_filename <- as.character(args[2])
# prom_or_enh_matrix_filename <- as.character(args[3])
# outname <- as.character(args[4])

setwd("/Users/naim/Desktop/plotting/per_locus/")
permutation_filename <- "region1_chr1_chisq_permuted"
observed_assoc_filename <- "region1_chr1_observed"
prom_or_enh_matrix_filename <- "region1_chr1_prom_or_enh_matrix.txt"
outname <- "region1_chr1_contrast_test_results2.txt"


library(data.table)

# Inputs:
#perm_assoc <- fread("region1_chr1_chisq_permuted", header=F, stringsAsFactors=T)
#obs_assoc <- fread("region1_chr1_observed", header=T, stringsAsFactors=F)
#prom_or_enh_matrix <- fread("region1_chr1_prom_or_enh_matrix.txt", header=T, stringsAsFactors=F)

print(paste("Permutation file:", permutation_filename))
print(paste("Observed association file:", observed_assoc_filename))
print(paste("Promoter or enhancer matrix file:", prom_or_enh_matrix_filename))
print("Reading data")
perm_assoc <- fread(permutation_filename, header=F, stringsAsFactors=T)
obs_assoc <- fread(observed_assoc_filename, header=T, stringsAsFactors=F)
prom_or_enh_matrix <- fread(prom_or_enh_matrix_filename, header=T, stringsAsFactors=F)


prom_enh_matrix <- as.matrix(prom_or_enh_matrix[,7:133,with=F])
obs_chi_matrix <- as.matrix(obs_assoc$Z^2)
perm_chi_matrix <- as.matrix(perm_assoc[,7:10006,with=F])
print(paste("Number of SNPs:", nrow(obs_chi_matrix)))
missings <- which(is.na(obs_chi_matrix[,1]) | is.na(perm_chi_matrix[,1]) | is.na(perm_chi_matrix[,2]) | is.na(perm_chi_matrix[,3]))
missings_dets <- perm_assoc[missings, 1:6,with=F]
print(paste("SNPs with missing values:", length(missings)))
print(paste("Details of SNPs being removed due to missing values:"))
print(missings_dets)
print(paste("Number of SNPs after removal of missing values:", nrow(obs_assoc) - length(missings)))
print("")
lowMAFindex <- which(perm_assoc[,6,with=F] < 0.01)
print(paste("Number of SNPs with MAF<1%:", length(lowMAFindex)))
print("Details of SNPs being removed due to low MAF")
lowMAFdets <- perm_assoc[lowMAFindex, 1:6, with=F]
print(lowMAFdets)
snps_to_remove <- unique(c(missings,lowMAFindex))
print(paste("Number of SNPs after removal of low MAF SNPs:", nrow(obs_assoc) - length(snps_to_remove)))

if(length(snps_to_remove)>0) {
  perm_assoc <- perm_assoc[-snps_to_remove,]
  obs_assoc <- obs_assoc[-snps_to_remove,]
  prom_or_enh_matrix <- prom_or_enh_matrix[-snps_to_remove,]
  prom_enh_matrix <- prom_enh_matrix[-snps_to_remove,]
  obs_chi_matrix <- obs_chi_matrix[-snps_to_remove,]
  perm_chi_matrix <- perm_chi_matrix[-snps_to_remove,]
}


# Outputs:
print("Matrix multiplication step:") 
print(paste("Promoter OR enhancer matrix dimension:", nrow(prom_enh_matrix), "rows (SNPs)", ncol(prom_enh_matrix),"columns (tissues)"))
print(paste("Permutations' matrix dimension:", nrow(perm_chi_matrix), "rows (SNPs)", ncol(perm_chi_matrix), "columns (permutations)"))
Chisum_obs_per_tissue <- t(prom_enh_matrix) %*% obs_chi_matrix
Chisum_perm_per_tissue <- t(prom_enh_matrix) %*% perm_chi_matrix

print("Now flipping 0s and 1s in promoter OR enhancer matrix for contrast")
flip <- function(v) {
  return(ifelse(v==0,1,0))
}
flipped_prom_enh_matrix <- apply(prom_enh_matrix, 2, flip)

flipped_chisum_obs_per_tissue <- t(flipped_prom_enh_matrix) %*% obs_chi_matrix
flipped_chisum_perm_per_tissue <- t(flipped_prom_enh_matrix) %*% perm_chi_matrix

contrast_test_per_tissue <- data.frame(tissue=colnames(prom_or_enh_matrix)[7:ncol(prom_or_enh_matrix)], contrast_stat=NA, empirical_p=NA, num_prom_enh_snps=NA, enh_prom_snps=NA)
contrast_test <- Chisum_obs_per_tissue/colSums(prom_enh_matrix) - flipped_chisum_obs_per_tissue/colSums(flipped_prom_enh_matrix)
contrast_test_null_dist <- Chisum_perm_per_tissue/colSums(prom_enh_matrix) - flipped_chisum_perm_per_tissue/colSums(flipped_prom_enh_matrix)
print("Calculating empirical P-values")
calculate_empP <- function(v) {
  obs=v[1]
  perm=v[2:length(v)]
  return(sum(perm>=obs) / (length(perm)+1))
}
tempMAT <- cbind(contrast_test, contrast_test_null_dist)
empP <- apply(tempMAT,1,calculate_empP)

num_enh_snps <- matrix(nrow=nrow(contrast_test_per_tissue),ncol=1)
enh_prom_snps_per_tissue <- matrix(nrow=nrow(contrast_test_per_tissue),ncol=1)
for(i in 1:nrow(contrast_test_per_tissue)) {
  enh_snps=""
  indexes <- which(prom_enh_matrix[,i]==1)
  num_enh_snps[i,1] <- length(indexes)
  enh_snps <- paste(as.matrix(perm_assoc[indexes,2,with=F]),collapse=",")
  enh_prom_snps_per_tissue[i,1] <- enh_snps
}

contrast_test_per_tissue$contrast_stat <- contrast_test
contrast_test_per_tissue$empirical_p <- empP
contrast_test_per_tissue$num_prom_enh_snps <- num_enh_snps
contrast_test_per_tissue$enh_prom_snps <- enh_prom_snps_per_tissue

print(paste("Writing results to",outname))
write.table(contrast_test_per_tissue, outname, quote=F, row.names=F, col.names=T)
print("Done")
print("") 
