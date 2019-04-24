#!/usr/bin/R
# Takes the autosomal and sex chromosome PLINK-calculated heterozyosity values and plots them
# Syntax: Rscript heterozygosity.R <plink_autosomes.het> <plink_sex_chr.het> <famfile.fam> <prefix (optional)>
# Outputs: hetAuto_vs_hetSex.pdf

args=(commandArgs(TRUE))

autosomes_het_file <- as.character(args[1])
sex_het_file <- as.character(args[2])
famfile <- as.character(args[3])
prefix <- as.character(args[4])
if(is.na(prefix)) {
  prefix <- ""
} else {
  prefix <- paste0(prefix,"_")
}

autosomes <- read.table(autosomes_het_file, header=T, stringsAsFactors=F)
XChr <- read.table(sex_het_file, header=T, stringsAsFactors=F)
sex_codes <- read.table(famfile, col.names=c("FID", "IID", "Parent1", "Parent2", "Sex", "Phenotype"))

het_autosome <- 1 - (autosomes$O.HOM. / autosomes$N.NM.)
het_XChr <- 1 - (XChr$O.HOM. / XChr$N.NM.)
autosomesSex <- merge(autosomes, sex_codes) # merge will match FID and IID columns as they have same name
XChrSex <- merge(XChr, sex_codes)

autosomeMales <- subset(autosomesSex, Sex %in% 1) 
autosomeFemales <- subset(autosomesSex, Sex %in% 2)
autosomeNoSex <- subset(autosomesSex, Sex %in% 0)
XMales <- subset(XChrSex, Sex %in% 1)
XFemales <- subset(XChrSex, Sex %in% 2)
XNoSex <- subset(XChrSex, !(Sex %in% 1 | Sex %in% 2))

HetAutoM <- data.frame(Het=1-(autosomeMales$O.HOM./autosomeMales$N.NM.),
                       FID=autosomeMales$FID, IID=autosomeMales$IID)
HetAutoF <- data.frame(Het=1-(autosomeFemales$O.HOM./autosomeFemales$N.NM.),
                       FID=autosomeFemales$FID, IID=autosomeFemales$IID)
HetSexM <- data.frame(Het=1-(XMales$O.HOM./XMales$N.NM.),
                     FID=XMales$FID, IID=XMales$IID)
HetSexF <- data.frame(Het=1-(XFemales$O.HOM./XFemales$N.NM.),
                     FID=XFemales$FID, IID=XFemales$IID)
HetAutoNoSex <- data.frame(Het=1-(autosomeNoSex$O.HOM./autosomeNoSex$N.NM.),
                          FID=autosomeNoSex$FID, IID=autosomeNoSex$IID)
HetSexNoSex <- data.frame(Het=1-(XNoSex$O.HOM./XNoSex$N.NM.),
                          FID=XNoSex$FID, IID=XNoSex)

# plot
xmin <- min(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)-min(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)/10
xmax <- max(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)+max(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)/10
ymin <- min(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)-min(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)/10
ymax <- max(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)+max(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)/10
xrange <- c(xmin,xmax)
yrange <- c(ymin,ymax)

pdf(paste0(prefix, "hetAuto_vs_hetSex.pdf"))
plot(xrange,yrange,xlab="Autosomal chromosome heterozygosity",
     ylab="Sex chromosome heterozygosity",type="n")
points(HetAutoM$Het, HetSexM$Het, col="blue")
points(HetAutoF$Het, HetSexF$Het, col="red")
points(HetAutoNoSex$Het, HetSexNoSex$Het, col="black", bg="black", pch=21)
legend("topleft", legend=c("Males", "Females", "Unknown"), col=c("blue", "red","black"), pch=c(21,21,19))
dev.off()

#boxplot(HetSexM$Het, main="Males Sex Chromosome Heterozygosity") # males should all be around zero
#boxplot(subset(HetSexM$Het, HetSexM$Het<0.03), main="Males Sex Chromosome Heterozygosity\nAfter Outlier Removal") # males should all be around zero

#boxplot(HetSexF$Het) # females should be above zero

# To create interactive plotly plot:
#males = data.frame(autoHet=HetAutoM$Het, sexHet=HetSexM$Het, sex="males", IID=HetAutoM$IID)
#females = data.frame(autoHet=HetAutoF$Het, sexHet=HetSexF$Het, sex="females", IID=HetAutoF$IID)
#unknowns = data.frame(autoHet=HetAutoNoSex$Het, sexHet=HetSexNoSex$Het, sex="unknown", IID=HetAutoNoSex$IID)
#data = rbind(males, females, unknowns)

#p <- plot_ly(data=data, x = ~autoHet, y = ~sexHet, color = ~sex, text = ~IID, type = 'scatter', mode = 'markers')
#chart_link = api_create(p, filename="BIOJUME-1903-Heterozygosity")
#chart_link


