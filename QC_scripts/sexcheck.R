#!/usr/bin/R
# Plots F statistic vs gender and sex obtained from a .sexcheck file from PLINK's --check-sex output
# Syntax: Rscript sexcheck.R <filename.sexcheck> [prefix]
# Outputs: (prefix)-sexcheck.pdf and png and a table of mismatches (prefix)-sexcheck.txt

args=(commandArgs(TRUE))

library(data.table)
library(ggplot2)


sexcheckfilename <- as.character(args[1])
prefix <- as.character(args[2])
if(is.na(prefix)) {
  prefix <- ""
} else {
  prefix <- paste0(prefix,"-")
}
sexcheck <- fread(sexcheckfilename, header=TRUE, stringsAsFactors=F)
sexcheck <- na.omit(sexcheck) # likely inconsequential
sexcheck$PEDSEX <- factor(sexcheck$PEDSEX, levels=0:2, labels=c("Unknown","Male","Female"))
sexcheck$SNPSEX <- factor(sexcheck$SNPSEX, levels=0:2, labels=c("Unknown","Male","Female"))
colnames(sexcheck)[which(colnames(sexcheck) %in% "PEDSEX")] <- "Gender"
colnames(sexcheck)[which(colnames(sexcheck) %in% "SNPSEX")] <- "Imputed_Sex"
#colnames(sexcheck)[which(colnames(sexcheck) %in% "F")] <- "chrXInbreeding coefficient (F)"

maleGender <- subset(sexcheck, sexcheck$Gender %in% "Male")
femaleGender <- subset(sexcheck, sexcheck$Gender %in% "Female")

pdf(paste0(prefix,"sexcheck.pdf"))
p <- ggplot(sexcheck, aes(x=Gender, y=F, color=STATUS, fill=STATUS))
#p + geom_boxplot()
#p + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6)
p + geom_jitter(shape=16, position=position_jitter(0.2), size=3) +
    geom_text(aes(label=ifelse(STATUS=="PROBLEM",as.character(IID),'')), position = position_jitter(0.2), hjust=0,vjust=1)
dev.off()
png(paste0(prefix,"sexcheck.png"),type="cairo")
p <- ggplot(sexcheck, aes(x=Gender, y=F, color=STATUS, fill=STATUS))
p + geom_jitter(shape=16, position=position_jitter(0.2), size=3) +
  geom_text(aes(label=ifelse(STATUS=="PROBLEM",as.character(IID),'')), position = position_jitter(0.2), hjust=0,vjust=1)
dev.off()

problematic <- subset(sexcheck, sexcheck$STATUS %in% "PROBLEM")
write.table(problematic, paste0(prefix,"sexcheck.txt"), row.names = F, col.names = T, quote = F, sep="\t")

