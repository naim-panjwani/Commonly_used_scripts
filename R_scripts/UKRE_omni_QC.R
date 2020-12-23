setwd("/Users/naim panjwani/Documents/Strug/UKRE_omni_QC_Illumina_calls")

imiss <- read.table("01_missing.imiss", header=TRUE)
lmiss <- read.table("01_missing.lmiss", header=TRUE)

Nsnps <- dim(lmiss)[1]
#pdf("imiss.pdf")
plot(imiss$F_MISS*100, xlab="Individual", ylab="Individual Missingness (%)")
#dev.off()
#plot of missing rate vs individual


high_imiss <- subset(imiss, F_MISS > 0.15)
IndHighMiss <- data.frame(high_imiss$FID, high_imiss$IID)
write.table(IndHighMiss, file="IndHighMiss.txt", quote=F, row.names=F, col.names=F)
#writing FID/IID combination of individuals with high SNP missing rate to file
#can use plink to exclude these individuals from analysis 

#setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3/imputation/trial2/")
#imiss <- read.table("Step5_missing.imiss", header=TRUE)
#lmiss <- read.table("Step5_missing.lmiss", header=TRUE)

plot((1-lmiss$F_MISS)*100, xlab="SNP", ylab="SNP call rate (%)")
Nindividuals <- dim(imiss)[1]
#lmiss2 <- subset(lmiss, lmiss$F_MISS > .02)
#complete_SNPs <- subset(lmiss, lmiss$F_MISS != 0)
#fully_missing_SNPs <- subset(lmiss, lmiss$F_MISS == 0)
lmiss10 <- subset(lmiss, lmiss$F_MISS > 0.10)
SNPmissList10 <- lmiss10$SNP
#SNPmissList2 <- lmiss2$SNP

write.table(SNPmissList10, file="SNPmissList10.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(SNPmissList2, file="SNPmissList2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#writing to file all SNPs that have a rate of > 0.10 missing across all 64 individuals
#could use this file to exclude these snps in plink



##### MAF
frqFile <- read.table("02_freq.frq", header=TRUE)
MAFsnps <- subset(frqFile, frqFile$MAF < 0.02)
library(lattice)

#pdf("MAFhist.pdf")
histogram(~frqFile$MAF,main="Minor Allele Frequency Histogram",
          xlab="Minor Allele Frequency")
#dev.off()
hist(frqFile$MAF,main="Minor Allele Frequency Histogram",
     xlab="MAF",ylab="Percent",freq=FALSE)

hist(frqFile$MAF,xlim=c(0,0.1), 
     main="MAF zoom-in",xlab="MAF",ylab="Percent",freq=FALSE)

# SNP file containing MAF < 2% (subset of ALL SNPs)
# Want to remove these using PLINK to analyze common variants
write.table(MAFsnps$SNP, file="MAFsnps.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)




#========================================================================================================
## Heterozygosity
## Sex check

autosomes <- read.table("09_UKREcnv_autosome_het.het", header=T)
XChr <- read.table("11_UKREcnv_pruned_chrX_het.het", header=T)
het_autosome <- 1 - (autosomes$O.HOM. / autosomes$N.NM.)
het_XChr <- 1 - (XChr$O.HOM. / XChr$N.NM.)

sex_codes <- read.table("10_UKREcnv_chr23to1_binary.fam", 
                        col.names=c("FID", "IID", "Parent1", "Parent2", "Sex", "Phenotype"))
autosomesSex <- merge(autosomes, sex_codes) # merge will match FID and IID columns as they have same name
XChrSex <- merge(XChr, sex_codes)

autosomeMales <- subset(autosomesSex, Sex %in% 1) 
autosomeFemales <- subset(autosomesSex, Sex %in% 2)
autosomeNoSex <- subset(autosomesSex, Sex %in% "Unknown")
XMales <- subset(XChrSex, Sex %in% 1)
XFemales <- subset(XChrSex, Sex %in% 2)
XNoSex <- subset(XChrSex, Sex %in% "Unknown")

HetAutoM <- data.frame(Het=1-(autosomeMales$O.HOM./autosomeMales$N.NM.),
                       FID=autosomeMales$FID, IID=autosomeMales$IID)
HetAutoF <- data.frame(Het=1-(autosomeFemales$O.HOM./autosomeFemales$N.NM.),
                       FID=autosomeFemales$FID, IID=autosomeFemales$IID)
HetSexM <- data.frame(Het=1-(XMales$O.HOM./XMales$N.NM.),
                       FID=XMales$FID, IID=XMales$IID)
HetSexF <- data.frame(Het=1-(XFemales$O.HOM./XFemales$N.NM.),
                     FID=XFemales$FID, IID=XFemales$IID)
HetAutoNoSex <- data.frame(Het=1-(autosomeNoSex$O.HOM./autosomeNoSex$N.NM.),
                          FID=autosomeNoSex$FID, IID=autosomeNoSex)
HetSexNoSex <- data.frame(Het=1-(XNoSex$O.HOM./XNoSex$N.NM.),
                         FID=XNoSex$FID, IID=XNoSex)

#plot
xmin <- min(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)-min(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)/10
xmax <- max(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)+max(HetAutoM$Het, HetAutoF$Het, HetAutoNoSex$Het)/10
ymin <- min(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)-min(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)/10
ymax <- max(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)+max(HetSexM$Het, HetSexF$Het, HetSexNoSex$Het)/10
xrange <- c(xmin,xmax)
yrange <- c(ymin,ymax)
plot(xrange,yrange,xlab="Autosomal chromosome heterozygosity",
     ylab="Sex chromosome heterozygosity",type="n")
points(HetAutoM$Het, HetSexM$Het, col="blue")
points(HetAutoF$Het, HetSexF$Het, col="red")
points(HetAutoNoSex$Het, HetSexNoSex$Het, col="black", bg="black", pch=21)
# each individual point represents one individual in dataset

plot(HetSexM$Het) # males should all be around zero
plot(HetSexF$Het) # females should be above zero

misclassified_males <- HetSexF[HetSexF$Het<0.02,]
write.csv(misclassified_males, "temp.csv", quote=F, row.names=F)


##### Plot of heterozygosity of autosomal chromosomes ############
##### (ie. due to possible sample contamination or actual inbreeding) ####
HetAuto <- rbind(HetAutoF,HetAutoM)
boxplot(HetAuto$Het, 
        main="Box Plot of Heterozygosity\nfor Autosomal Chromosomes")
iqr <- quantile(HetAuto$Het,0.75)-quantile(HetAuto$Het,0.25)
upperWhisk <- quantile(HetAuto$Het,0.75)+1.5*iqr
lowerWhisk <- quantile(HetAuto$Het,0.25)-1.5*iqr
#hetOutliers <- HetAuto[HetAuto$Het>upperWhisk | HetAuto$Het<lowerWhisk,]
#write.table(hetOutliers,file="02_hetOutliers2.het",quote=FALSE, row.names=FALSE)

(obvious_outliers <- HetAuto[HetAuto$Het<.28,]) #cut-off of 0.28 according to Jiafen 
# IID DOC is around 0.27
write.table(obvious_outliers[,2:3], file="hetOutliers.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# subset(HetAuto, grepl("7066",HetAuto$IID))
# The other 2 individuals in the same family 7066 have heterozygosity of 0.349-0.351

#(upper_outliers <- subset(HetAuto, HetAuto$Het>.4))
(lower_outliers <- subset(HetAuto, HetAuto$Het< 0.28))
for (i in 1:length(lower_outliers$Het)) points(lower_outliers$Het[i],col="red", pch="*",cex=1.5)
#for (i in 1:length(upper_outliers$Het)) points(upper_outliers$Het[i],col="red", pch="*",cex=1.5)









#========================================================================================================
# #Sex Check
# #after running check-sex option on plink, loading the .sexcheck output file
# #PLINK option declares males and females according to heterozygosity rate (F value)

sexCheck <- read.table("12_UKREcnv_checksex.sexcheck", header=TRUE)
sexErr <- with(sexCheck, which((PEDSEX != SNPSEX) | (F<0.8 & PEDSEX==1) | (F>0.2 & PEDSEX == 2)))
indWithSexProblem <- sexCheck[sexErr,]
indWithSexProblem2 <- subset(sexCheck, sexCheck$STATUS %in% "PROBLEM")

gender_mismatch <- subset(sexCheck, with(sexCheck, PEDSEX!=SNPSEX & SNPSEX!=0))
gender_mismatch_index <- which(with(sexCheck, PEDSEX!=SNPSEX & SNPSEX!=0))

F_values_to_plot <- sexCheck[which(with(sexCheck, PEDSEX!=1 & F<0.95)),]$F

boxplot(F_values_to_plot)
points(indWithSexProblem$F, col="red", pch="*",cex=1.5)
sapply(indWithSexProblem$F, points, col="red", pch="*",cex=1.5)


