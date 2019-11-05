#!/usr/bin/R

library(argparser, quietly=TRUE)
library(data.table)

p <- arg_parser("Given a gene and tissue, plot the expression values")
p <- add_argument(p, "expression_file", help="The expression file in TPM")
p <- add_argument(p, "gene", help="The gene name in ENSG or regular name")
p <- add_argument(p, "tissueName", help="The name of the tissue for the figure")
argv <- parse_args(p)

#expression_file="../GTEx_eQTL_Pipeline/output/02-RNAseq_Ileum_subset/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_Ileum.gz"
#gene="ENSG00000075673.7"
#outfilename="Ileum_tpm_ATP12A"

print(paste0("expression_file: ", argv$expression_file))
print(paste0("gene: ", argv$gene))
print(paste0("tissueName: ", argv$tissueName))

geneExpression <- fread(paste0("gunzip -c ", argv$expression_file), header=T, stringsAsFactors=F)
if(nrow(subset(geneExpression, geneExpression$Description %in% argv$gene)) > 0) {
  expValues <- subset(geneExpression, geneExpression$Description %in% argv$gene)
} else if (nrow(subset(geneExpression, geneExpression$Name%in% argv$gene)) > 0) {
  expValues <- subset(geneExpression, geneExpression$Name %in% argv$gene)
} else {
   stop("Gene not found")
}

print(paste0("Saving plot in: ", as.character(expValues[1,2,with=F]), "_", argv$tissueName, ".pdf"))

aboveThresh=sum(as.numeric(expValues[,3:ncol(expValues),with=F]) > 0.1)
pct = aboveThresh / (ncol(expValues)-2)
pdf(paste0(as.character(expValues[1,2,with=F]), "_", argv$tissueName, ".pdf"))
boxplot(as.numeric(expValues[,3:ncol(expValues),with=F]), ylab="Expression TPM",
        main=paste0(as.character(expValues[1,2,with=F])," ",argv$tissueName, " V7\n", "N=",(ncol(expValues)-2),
                    "\nAbove threshold: ", aboveThresh, "/", (ncol(expValues)-2), "=", pct) ,
        outpch=NA)
stripchart(as.numeric(expValues[,3:ncol(expValues),with=F]), vertical=TRUE, method="jitter", add=TRUE, pch=21, col="maroon", bg="bisque")
abline(a=0.1,b=0, col="red")
dev.off()

