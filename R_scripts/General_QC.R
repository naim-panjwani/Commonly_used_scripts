##### Missing rates of SNPs across individuals
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v3")

imiss <- read.table("Step2_missingness.imiss", header=TRUE)
lmiss <- read.table("Step2_missingness.lmiss", header=TRUE)

Nsnps <- dim(lmiss)[1]
#pdf("imiss.pdf")
plot(imiss$F_MISS*100, xlab="Individual", ylab="Individual Missingness (%)")
#dev.off()
#plot of missing rate vs individual


high_imiss <- subset(imiss, F_MISS > 0.1)
IndHighMiss <- data.frame(high_imiss$FID, high_imiss$IID)
write.table(IndHighMiss, file="IndHighMiss.txt", quote=F, row.names=F, col.names=F)
#writing FID/IID combination of individuals with high SNP missing rate to file
#can use plink to exclude these individuals from analysis 

#setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3/imputation/trial2/")
#imiss <- read.table("Step5_missing.imiss", header=TRUE)
#lmiss <- read.table("Step5_missing.lmiss", header=TRUE)

plot((1-lmiss$F_MISS)*100, xlab="SNP", ylab="SNP call rate (%)")
hist((1-lmiss$F_MISS)*100, main="SNP Call Rate", xlab="SNP Call Rate")

Nindividuals <- dim(imiss)[1]
lmiss2 <- subset(lmiss, lmiss$F_MISS > .02)
#complete_SNPs <- subset(lmiss, lmiss$F_MISS != 0)
#fully_missing_SNPs <- subset(lmiss, lmiss$F_MISS == 0)
lmiss10 <- subset(lmiss, lmiss$F_MISS > 0.10)
SNPmissList10 <- lmiss10$SNP
SNPmissList2 <- lmiss2$SNP

write.table(SNPmissList10, file="SNPmissList10.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(SNPmissList2, file="SNPmissList2.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
#writing to file all SNPs that have a rate of > 0.10 missing across all 64 individuals
#could use this file to exclude these snps in plink



##### MAF
frqFile <- read.table("Step6_maf.frq", header=TRUE)
MAFsnps <- subset(frqFile, frqFile$MAF < 0.02)
library(lattice)

#pdf("MAFhist.pdf")
histogram(~frqFile$MAF,main="Minor Allele Frequency Histogram",
          xlab="Minor Allele Frequency")
#dev.off()

hist(frqFile$MAF, main="MAF", xlab="MAF",ylab="Percent",freq=FALSE)

hist(frqFile$MAF,xlim=c(0,0.1), ylim=c(0,25),
     main="MAF zoom-in",xlab="MAF",ylab="Percent",freq=FALSE)

# SNP file containing MAF < 2% (subset of ALL SNPs)
# Want to remove these using PLINK to analyze common variants
write.table(MAFsnps$SNP, file="MAFsnps.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)



###### Heterozygous haploid SNPs; list in txt file to exclude using PLINK
# hh <- read.table("Step5_exclude_high_SNP_missing_rate.hh", header=FALSE)
# hhSNPs <- hh$V3
# write.table(hhSNPs, file="hh.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Can do within Unix


###### Heterozygosity/inbreeding check

#Heterozygosity--converting chromosome 23 to chromosome 1

# Recode binary PLINK files back to ped/map 
# then change the chromosome column to an autosomal chromosome
# in the map file; also rename the ped file using the mv command
##Chr23 <- read.table("Step7_chr23_pruned_data.map")
##Chr23$V1 <- 1
##write.table(Chr23, "Step7_chr23tochr1.map", quote=FALSE, row.names=FALSE, col.names=FALSE)
# renamed ped file to Step7_chr23tochr1.ped as well

# Now we can run the --het option in PLINK on X chromosome and verify gender
# Should also run --het option on unchanged data 
# (autosomal chromosomes) for inbreeding check later

autosomes <- read.table("09_2_autosomes.het", header=T)
XChr <- read.table("09_2_sexChr.het", header=T)
het_autosome <- 1 - (autosomes$O.HOM. / autosomes$N.NM.)
het_XChr <- 1 - (XChr$O.HOM. / XChr$N.NM.)

sex_codes <- read.table("09_pruned.fam", 
                        col.names=c("FID", "IID", "Parent1", "Parent2", "Sex", "Phenotype"))
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
                          FID=XNoSex$FID, IID=XNoSex$IID)

# plot
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
#each individual point represents one individual in dataset

plot(HetSexM$Het) # males should all be around zero
plot(HetSexF$Het) # females should be above zero

misclassified_males <- HetSexF[HetSexF$Het<0.02,]
write.csv(misclassified_males, "misclassified_males.csv", quote=F, row.names=F)


##### Plot of heterozygosity of autosomal chromosomes ############
##### (ie. due to possible sample contamination or actual inbreeding) ####
HetAuto <- rbind(HetAutoF,HetAutoM, HetAutoNoSex)
boxplot(HetAuto$Het, 
        main="Box Plot of Heterozygosity\nfor Autosomal Chromosomes")
iqr <- quantile(HetAuto$Het,0.75)-quantile(HetAuto$Het,0.25)
upperWhisk <- quantile(HetAuto$Het,0.75)+1.5*iqr
lowerWhisk <- quantile(HetAuto$Het,0.25)-1.5*iqr
#hetOutliers <- HetAuto[HetAuto$Het>upperWhisk | HetAuto$Het<lowerWhisk,]
#write.table(hetOutliers,file="02_hetOutliers2.het",quote=FALSE, row.names=FALSE)

(obvious_outliers <- HetAuto[HetAuto$Het<.28,]) #cut-off of 0.28 according to Jiafen 
# IID ? 
write.table(obvious_outliers[,2:3], file="hetOutliers.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

# subset(HetAuto, grepl("7066",HetAuto$IID))
# The other 2 individuals in the same family 7066 have heterozygosity of 0.349-0.351

(upper_outliers <- subset(HetAuto, HetAuto$Het>.4))
(lower_outliers <- subset(HetAuto, HetAuto$Het< 0.28))
#for (i in 1:length(lower_outliers$Het)) points(lower_outliers$Het[i],col="red", pch="*",cex=1.5)
#for (i in 1:length(upper_outliers$Het)) points(upper_outliers$Het[i],col="red", pch="*",cex=1.5)



#========================================================================================================
# #Sex Check
# #after running check-sex option on plink, loading the .sexcheck output file
# #PLINK option declares males and females according to heterozygosity rate (F value)
# 
# sexCheck <- read.table("Step10check_sex.sexcheck", header=TRUE)
# sexErr <- with(sexCheck, which((PEDSEX != SNPSEX) | (F<0.8 & PEDSEX==1) | (F>0.2 & PEDSEX == 2)))
# indWithSexProblem <- sexCheck[sexErr,]
# indWithSexProblem2 <- subset(sexCheck, sexCheck$STATUS %in% "PROBLEM")
# # 7062_202 is being marked as PROBLEM b/c F>0.2 (F=0.2207); will ignore for now
# # as it is only slightly above 0.2
# females <- subset(sexCheck, (SNPSEX %in% 2) | (PEDSEX %in% 2))
# boxplot(females$F)



#========================================================================================================
############### HWE analysis ###################################
hwe <- read.table("Step11hardy_on_all_subjects.hwe", header=TRUE)
#unaffected <- subset(hwe, TEST %in% "UNAFF")
#expectedP <- runif(n=length(unaffected$P),0,1)  #expected p-values are uniformly-distributed
#orderedP <- unaffected[order(unaffected$P),]
#qqplot(-log10(expectedP), -log10(orderedP$P), main="QQ Plot of Observed vs Expected HWE p-values", xlab="-log10(Expected p-value) quantile", ylab="-log10(Observed HWE p-value) quantile")
#abline(0,1,col="red")
#abline(-log10(5*10^(-7)),0,col="red")

orderedHWE <- hwe[order(hwe$P),]

filter1 <- 5*10^(-7)
#filter1hwe <- orderedP[orderedP$P>filter1,]
#exclude1hwe <- orderedP[orderedP$P<filter1,]
filter1hweTEMP <- orderedHWE[orderedHWE$P>filter1,]
exclude1hweTEMP <- orderedHWE[orderedHWE$P<filter1,]

filter2 <- 10^(-4)
#filter2hwe <- orderedP[orderedP$P>filter2,]
#exclude2hwe <- orderedP[orderedP$P<filter2,]
filter2hweTEMP <- orderedHWE[orderedHWE$P>filter2,]
exclude2hweTEMP <- orderedHWE[orderedHWE$P<filter2,]

filter3 <- 10^(-3)
filter3hwe <- orderedP[orderedP$P>filter3,]
exclude3hwe <- orderedP[orderedP$P<filter3,]

filter4 <- 10^(-2)
filter4hwe <- orderedP[orderedP$P>filter4,]
exclude4hwe <- orderedP[orderedP$P<filter4,]

# write.table(filter1hwe, "09_hweFilter1SNPs.hwe",quote=FALSE,row.names=FALSE)
# write.table(filter2hwe, "09_hweFilter2SNPs.hwe",quote=FALSE,row.names=FALSE)
# write.table(filter3hwe, "09_hweFilter3SNPs.hwe",quote=FALSE,row.names=FALSE)
# write.table(filter4hwe, "09_hweFilter4SNPs.hwe",quote=FALSE,row.names=FALSE)
# write.table(exclude1hwe, "09_hweExclude1SNPs.hwe",quote=FALSE,row.names=FALSE)
write.table(exclude2hweTEMP$SNP, "HWE_on_all_subjects.txt",
            quote=FALSE,row.names=FALSE, col.names=FALSE)
# write.table(exclude3hwe, "09_hweExclude3SNPs.hwe",quote=FALSE,row.names=FALSE)
# write.table(exclude4hwe, "09_hweExclude4SNPs.hwe",quote=FALSE,row.names=FALSE)




#========================================================================================================
# Cryptic relatedness/duplicates

# Read in genome file generated by PLINK after re-pruning the data
relatedness <- read.table("12_kinship.genome", header=TRUE)
library(rgl)
with(relatedness, plot3d(Z0,Z1,Z2))

### 1. Mark the different groups according to their (Z0,Z1,Z2) coordinates
#unrelated <- subset(relatedness, with(relatedness, Z0>0.3 & Z2<0.1)) 
#includes first cousins, half-sibs & double first cousins

#################### CHANGED CRITERIA FOR UNRELATED / 1ST COUSINS / HALF SIBS 
#################### USE CODE IN NEXT SECTION

# unrelated2 <- subset(relatedness, with(relatedness, Z0>0.9 & Z2<0.1))
#   #excludes first cousins, half-sibs and double first cousins
# 
# first_cousins <- subset(relatedness, with(relatedness, Z0>0.7 & Z0<=0.9 & Z2<0.1))
# halfsib <- subset(relatedness, with(relatedness, Z0>0.3 & Z0<=0.7 & Z2<0.1))
# sib_pairs <- subset(relatedness, with(relatedness, Z0>0.1 & Z1>0.3 & Z2>0.1))
# parent_offspring <- subset(relatedness, with(relatedness, Z0<0.1 & Z1>0.8 & Z2<0.2))
# twins_duplicates <- subset(relatedness, with(relatedness, Z2>0.9))
# 
# # ensure all pairs have been included:
# (all_included <- dim(relatedness)[1] == (dim(unrelated2)[1]+dim(first_cousins)[1]+
#                                            dim(sib_pairs)[1]+dim(halfsib)[1]+
#                                           dim(parent_offspring)[1]+dim(twins_duplicates)[1]) )
# # TRUE
# 
# with(relatedness, plot3d(Z0,Z1,Z2, type="n"))
# with(unrelated2, points3d(Z0,Z1,Z2,col="black",size=4))
# with(first_cousins, points3d(Z0,Z1,Z2,col="yellow",size=4))
# with(halfsib, points3d(Z0,Z1,Z2,col="magenta",size=4))
# with(sib_pairs, points3d(Z0,Z1,Z2,col="green",size=4))
# with(parent_offspring, points3d(Z0,Z1,Z2,col="blue",size=4))
# with(twins_duplicates, points3d(Z0,Z1,Z2,col="red",size=5))
# 
# ### 2. Verify family grouping is marked correctly in dataset
# # a) For different FIDs, identify if any relative pairs have been missed
# differing_FID <- subset(relatedness, relatedness$FID1 != relatedness$FID2)
# unmarked_relatives <- subset(differing_FID, with(differing_FID, Z0<0.3 & Z2>0.1))
# dim(unmarked_relatives)
# # 9 pairs that have been missed or marked in the incorrect family
# with(unmarked_relatives, points3d(Z0,Z1,Z2, col="red",size=6))
# print(unmarked_relatives[,c(1:5,7:10)] )
# # some are sibling pairs and others are twins/duplicates
# 
# #### Group the relative pairs into the same family (FIDs should be equal)
# # Load fam file and correct FIDs for relatives
# fam <- read.table("Step13_pruned.fam", header=FALSE)
# #unmarked_relatives_siblings <- subset(unmarked_relatives, unmarked_relatives$PI_HAT<0.9)
# #unmarked_relatives_duplicates <- subset(unmarked_relatives, unmarked_relatives$PI_HAT>0.9)
# 
# # Find the misclassified relatives (FID2 IID2)
# indices <- which(fam$V2 %in% unmarked_relatives$IID2)
# # replace FID in fam file with FID1 in unmarked_relatives
# for (i in 1:dim(fam[indices,])[1]) {
#   IID <- fam[indices[i],]$V2
#   IID2index <- which(unmarked_relatives$IID2 %in% IID)
#   fam[indices[i],]$V1 <- unmarked_relatives[IID2index,"FID1"]
# }
# # Verify correct replacement:
# fam[indices,]
# 
# # Create new updated fam file
# write.table(fam, "Step13_pruned_updatedFIDs.fam",
#             quote=FALSE,row.names=FALSE, col.names=FALSE)
# 
# length(unique(sib_pairs$IID1))
# length(unique(parent_offspring$IID1))
# length(unique(twins_duplicates$IID1))
# 
# 
# ### Next verify that IBD for members of the same family make sense
# familyIBD <- read.table("Step15_rel_check.genome", header=TRUE)
# with(familyIBD, plot3d(Z0,Z1,Z2, col="magenta",size=3))
# 
# 
# ## 3. For pairs with the SAME FID
# fam_unrelated <- NULL
# for (i in unique(relatedness$FID1)) {
#   family <- subset(relatedness, relatedness$FID1 == i & relatedness$FID2 == i)
#   unrelateds <- subset(family, with(family,  Z0>0.9 & Z2<0.1)) 
#     #cousins/half-sibs assumed unrelated
#   fam_unrelated <- rbind(fam_unrelated, unrelateds)
# }
# dim(fam_unrelated)
# length(unique(fam_unrelated$IID1))
# 
# # write.table(sib_pairs, "siblings.txt",
# #             quote=FALSE,row.names=FALSE, col.names=TRUE)
# # write.table(parent_offspring, "parent_offspring.txt",
# #             quote=FALSE,row.names=FALSE, col.names=TRUE)
# # write.table(first_cousins, "first_cousins.txt",
# #             quote=FALSE,row.names=FALSE, col.names=TRUE)
# # write.table(halfsib, "half_sibs.txt",
# #             quote=FALSE,row.names=FALSE, col.names=TRUE)
# # write.table(twins_duplicates, "duplicates.txt",
# #             quote=FALSE,row.names=FALSE, col.names=TRUE)

#===================================================================================
# Refining unrelated 2nd and 3rd degree relatives further
library(lattice)
library(rgl)
#relatedness <- read.table("03_hapmap_hg19_kinship.genome", header=TRUE)
relatednessCopy<-relatedness
with(relatednessCopy, plot3d(Z0,Z1,Z2))

relatedness<-subset(relatedness, relatedness$PI_HAT!=0)


print(densityplot(~relatedness$PI_HAT, width=0.04)) # overall plot overwhelmed by unrelated individuals
print(densityplot(~relatedness$PI_HAT, from=0.15, to=0.25,width=0.04)) # some 3rd-degree relatives (ie. first cousins) but hard to resolve
print(densityplot(~relatedness$PI_HAT, from=0.09, to=0.7,width=0.04))

# Let's consider the relatedness groups a little differently...

unrelated2 <- subset(relatedness, with(relatedness, Z2<0.23 & PI_HAT < 0.24))
with(relatedness, plot3d(Z0,Z1,Z2, type="n"))
with(unrelated2, points3d(Z0,Z1,Z2, col="black"))
unrelated2_pairs <- paste(as.character(unrelated2$IID1), as.character(unrelated2$IID2),sep=",")
# first cousins considered unrelated

#first_cousins <- subset(relatedness, with(relatedness, Z0>0.7 & Z0<=0.9 & Z2<0.1))
second_degree_relatives <- subset(relatedness, with(relatedness, Z0>0.3 & Z2<0.1 & 
                                                      PI_HAT >=0.24 & PI_HAT <0.4))
second_degree_pairs <- paste(as.character(second_degree_relatives$IID1), as.character(second_degree_relatives$IID2),sep=",")
with(second_degree_relatives, points3d(Z0,Z1,Z2, col="darkgoldenrod3"))
sib_pairs <- subset(relatedness, with(relatedness, Z0>0.1 & Z1>0.3 & Z2>=0.1))
sib_pairs_pairs <- paste(as.character(sib_pairs$IID1), as.character(sib_pairs$IID2),sep=",")
with(sib_pairs, points3d(Z0,Z1,Z2, col="blue"))
parent_offspring <- subset(relatedness, with(relatedness, Z0<0.1 & Z1>0.8 & Z2<0.2))
parent_offspring_pairs <- paste(as.character(parent_offspring$IID1), as.character(parent_offspring$IID2),sep=",")
with(parent_offspring, points3d(Z0,Z1,Z2, col="red"))
twins_duplicates <- subset(relatedness, with(relatedness, Z2>0.9))
duplicate_pairs <- paste(as.character(twins_duplicates$IID1), as.character(twins_duplicates$IID2),sep=",")
with(twins_duplicates, points3d(Z0,Z1,Z2, col="darkgreen"))

(all_included <- dim(relatedness)[1] == (dim(unrelated2)[1]+dim(second_degree_relatives)[1]+dim(sib_pairs)[1]+
                                           dim(parent_offspring)[1]+dim(twins_duplicates)[1]) )
# FALSE

# included <- rbind(unrelated2, second_degree_relatives, sib_pairs, parent_offspring, twins_duplicates)
# unique_pairs_included <-  paste(as.character(included$IID1),as.character(included$IID2),sep=",")
# all_unique_pairs <- paste(as.character(relatedness$IID1), as.character(relatedness$IID2),sep=",")
# who_is_missing_index <- which(!(all_unique_pairs %in%  unique_pairs_included))
# who_is_missing <- all_unique_pairs[who_is_missing_index]
# which(all_unique_pairs %in% who_is_missing)
# relatedness[which(all_unique_pairs %in% who_is_missing),]
#dups <- unique_pairs_included[duplicated(unique_pairs_included)]
#sib_pairs[which(sib_pairs_pairs %in% dups),]
#second_degree_relatives[which(second_degree_pairs %in% dups),]


write.table(second_degree_relatives, "second_degree_relatives.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(sib_pairs, "siblings.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(parent_offspring, "parent_offspring.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(twins_duplicates, "duplicates.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)
write.table(unrelated2, "unrelated.txt",
            quote=FALSE,row.names=FALSE, col.names=TRUE)

with(relatedness, plot3d(Z0,Z1,Z2, type="n"))
with(unrelated2, points3d(Z0,Z1,Z2,col="black",size=4))
with(second_degree_relatives, points3d(Z0,Z1,Z2,col="yellow",size=4))
with(sib_pairs, points3d(Z0,Z1,Z2,col="green",size=4))
with(parent_offspring, points3d(Z0,Z1,Z2,col="blue",size=4))
with(twins_duplicates, points3d(Z0,Z1,Z2,col="red",size=5))

#=============================================================== Kinship cont'd
# Remove one related pair from sample based on missingness
imiss <- read.table("02_missing.imiss", header=TRUE)


#------------------------------------------------------------------------------------------------
WorstPair <- function(related_pairs, imiss) {
  # PRE: takes in .genome format data.frame and .imiss missingness information
  # POST: returns vector of each pair in .genome file with highest missingness
  
  num_pairs <- dim(related_pairs)[1]
  to_be_removed <- NULL
  for (i in 1:num_pairs) { # for each pair
    miss1 <- subset(imiss, imiss$IID %in% related_pairs$IID1[i])$F_MISS
    miss2 <- subset(imiss, imiss$IID %in% related_pairs$IID2[i])$F_MISS
    if(miss1 > miss2) { 
      pair_to_be_removed <- cbind(as.character(related_pairs$FID1[i]),as.character(related_pairs$IID1[i])) 
    } else {
      pair_to_be_removed <- cbind(as.character(related_pairs$FID2[i]),as.character(related_pairs$IID2[i])) 
    }
    to_be_removed <- rbind(to_be_removed,pair_to_be_removed)
  }
  return(to_be_removed)
}
#------------------------------------------------------------------------------------------------
PairClassify <- function(related_pairs, imiss) {
  # PRE: takes in .genome format data.frame and .imiss missingness information
  # POST: returns vector of each pair in .genome file with highest missingness
  
  num_pairs <- dim(related_pairs)[1]
  to_be_removed <- NULL
  comp_pairs <- NULL
  for (i in 1:num_pairs) { # for each pair
    miss1 <- subset(imiss, imiss$IID %in% related_pairs$IID1[i])$F_MISS
    miss2 <- subset(imiss, imiss$IID %in% related_pairs$IID2[i])$F_MISS
    if(miss1 > miss2) { 
      pair_to_be_removed <- cbind(as.character(related_pairs$FID1[i]),as.character(related_pairs$IID1[i])) 
      complementary_pair <- cbind(as.character(related_pairs$FID2[i]),as.character(related_pairs$IID2[i])) 
    } else {
      pair_to_be_removed <- cbind(as.character(related_pairs$FID2[i]),as.character(related_pairs$IID2[i])) 
      complementary_pair <- cbind(as.character(related_pairs$FID1[i]),as.character(related_pairs$IID1[i])) 
    }
    to_be_removed <- rbind(to_be_removed,pair_to_be_removed)
    comp_pairs <- rbind(comp_pairs,complementary_pair)
  }
  return(list(to_be_removed, comp_pairs))
}
#------------------------------------------------------------------------------------------------
#### Bad way of removing pairs because the related pair may already
#### had been removed and then you are removing another one….
#
#pairs_to_remove <- rbind(WorstPair(second_degree_relatives,imiss),
#                     WorstPair(sib_pairs,imiss),
#                     WorstPair(parent_offspring,imiss),
#                     WorstPair(twins_duplicates,imiss))
#best_pair <- rbind(PairClassify(second_degree_relatives,imiss)[[2]],
#                   PairClassify(sib_pairs,imiss)[[2]],
#                   PairClassify(parent_offspring,imiss)[[2]],
#                   PairClassify(twins_duplicates,imiss)[[2]])
#

#write.table(pairs_to_remove, "pairs_to_remove.txt",
#            )

#------------------------------------------------------------------------------------------------
countDuplicates <- function(x, v) {
  k<-0
  for (i in 1:length(v)) {
    if(x==v[i]) {
      k<-k+1
    }
  }
  return(k)
}
#------------------------------------------------------------------------------------------------
OptimizeRelatedRemoval <- function(related_pairs, imiss) {
  to_be_removed <- NULL
  fids<-NULL
  j <- 1
  for (i in 1:dim(related_pairs)[1]) { # for each pair
    if(!(as.character(related_pairs$IID1[i]) %in% to_be_removed) & !(as.character(related_pairs$IID2[i]) %in% to_be_removed)) {
      current_FID <- as.character(related_pairs$FID1[i])
      if(current_FID != as.character(related_pairs$FID2[i])) stop("Relateds not in same family! FID1: ",related_pairs$FID1,"FID2: ",related_pairs$FID2)
      family_size <- countDuplicates(current_FID, as.character(related_pairs$FID1))
      if(countDuplicates(as.character(related_pairs$IID1[i]),c(as.character(related_pairs$IID1),as.character(related_pairs$IID2))) > 
           countDuplicates(as.character(related_pairs$IID2[i]), c(as.character(related_pairs$IID1),as.character(related_pairs$IID2)))) {
        to_be_removed[j] <- as.character(related_pairs$IID1[i])
        fids[j]<-as.character(related_pairs$FID1[i])
        j<-j+1
      } else if (countDuplicates(as.character(related_pairs$IID1[i]),c(as.character(related_pairs$IID1),as.character(related_pairs$IID2))) == 
                   countDuplicates(as.character(related_pairs$IID2[i]), c(as.character(related_pairs$IID1),as.character(related_pairs$IID2)))) {
        temp <- WorstPair(related_pairs[i,],imiss)
        to_be_removed[j] <- temp[1,2]
        fids[j] <- temp[1,1]
        j<-j+1
      } else {
        to_be_removed[j] <- as.character(related_pairs$IID2[i])
        fids[j]<-as.character(related_pairs$FID2[i])
        j<-j+1
      }
    }
  }
  return(as.data.frame(cbind(FID=as.character(fids),IID=as.character(to_be_removed))))
}


related <- rbind(second_degree_relatives, sib_pairs, parent_offspring, twins_duplicates)
related_individuals_to_remove <- WorstPair(related, imiss)
temp <- OptimizeRelatedRemoval(related, imiss)

write.table(related_individuals_to_remove, "related_individuals_to_remove.txt",
            quote=F, row.names=F, col.names=F)

# Now use PLINK to remove related individuals written to file






#============================================================================
# Stratify by ethnicity


#==================== IGNORE THIS CODE; TAKES TOO LONG ===================#
# 1. Identify common SNPs in Hapmap file
# hapmap <- read.table("hapmap_with_hg19_cor.bim", header=FALSE)
# UKRE <- read.table("Step9_hetOutliers_removed.bim", header=FALSE)
# hapmapSNPs <- hapmap$V2
# UKREsnps <- UKRE$V2
# rm(hapmap)
# rm(UKRE)
# 
# #install.packages("R.utils")
# #library(R.utils)
# commonSNPs <- NULL
# #total <- length(UKREsnps)
# #pb <- txtProgressBar(min = 0, max = total, style = 3) #create progress bar
# for (i in 1:length(UKREsnps)) {
#   #Sys.sleep(0.1) #update progress bar
#   #setTxtProgressBar(pb, i)
#   snp_is_common <- FALSE
#   j=1
#   #totalj <- length(hapmapSNPs)
#   #pbj <- txtProgressBar(min = 0, max = totalj, style = 3)
#   while (!snp_is_common & j != length(hapmapSNPs)+1) {
#     if (as.character(UKREsnps[i]) == as.character(hapmapSNPs[j])) {
#       commonSNPs <- c(commonSNPs,paste(as.character(UKREsnps[i]),sep="\n"))
#       snp_is_common <- TRUE
#     }
#     j <- j+1
#     if (j==length(hapmapSNPs) | snp_is_common) print(i)
#     #Sys.sleep(0.1)
#     #setTxtProgressBar(pbj,j)
#   }
#   #close(pbj)
# }
# #close(pb)
# 
# write.table(commonSNPs, "commonSNPs.txt",
#             quote=FALSE,row.names=FALSE, col.names=FALSE)


#====================================================================
# Stratify by ethnicity (better way)

# write the SNPs from RE dataset
# UKRE <- read.table("Step9_hetOutliers_removed.bim", header=FALSE)
# UKREsnps <- UKRE$V2
# write.table(UKREsnps, "UKREsnps.txt",
#             quote=FALSE,row.names=FALSE, col.names=FALSE)
# 
# # Use PLINK's extract command to extract the common SNPs from hapmap file
# # Write this list of common SNPs from hapmap file the extract these 
# #   from UKRE dataset to be left with just the common SNPs
# commonSNPs <- read.table("hapmap_common_snps.bim")$V2
# write.table(commonSNPs, "commonSNPs.txt",
#             quote=FALSE,row.names=FALSE, col.names=FALSE)

# Use PLINK to merge the hapmap and UKRE files
# Prune the merged file
# Remove related individuals

#====================================================================
# FUNCTIONS
#====================================================================
boxArea <- function(dataPoints) {
  # PRE: given 3 vectors of data points like 3 PCs in data.frame format
  # POST: returns the xrange, yrange and zrange of the 3D box as a vector 
  #     in the following format: c(xmin, xmax, ymin, ymax, zmin, zmax)
  xmin <- min(dataPoints[,1])
  xmax <- max(dataPoints[,1])
  ymin <- min(dataPoints[,2])
  ymax <- max(dataPoints[,2])
  zmin <- min(dataPoints[,3])
  zmax <- max(dataPoints[,3])
  return(c(xmin, xmax, ymin, ymax, zmin, zmax))
}

isInRange <- function(range, point) (point >= range[1] & point <= range[2])


isInBox <- function(boxPoints, samplePoints) {
  # PRE: give the 3 PCs from hapmap data for specific ethnic group in boxPoints and 
  #       the sample individuals to be tested for belonging in the
  #       cluster under samplePoints
  # POST: returns boolean vector on whether individual is within ethnic group
  inBox <- logical(length=dim(samplePoints)[1])
  box <- boxArea(boxPoints)
  xrange <- c(box[1], box[2])
  yrange <- c(box[3], box[4])
  zrange <- c(box[5], box[6])
  for (i in 1:dim(samplePoints)[1]) {
    if (isInRange(xrange, samplePoints[i,1]) & 
          isInRange(yrange, samplePoints[i,2]) &
          isInRange(zrange, samplePoints[i,3]))
    {
      inBox[i] <- TRUE
    } else {
      inBox[i] <- FALSE
    }
  }
  return(inBox)
}
#====================================================================
# END OF FUNCTIONS
#====================================================================

# Plot screeplot and individuals' eigenvectors

# read files
setwd("/Users/naim/Documents/Strug/UKRE/UKRE_QC/exomechip/UKRE_QC_v3")
#setwd("/Users/naim panjwani/Documents/Strug/UKRE_QC_v3")

eigenvalues <- read.table("44_dataset_hapmap_autosomalSNPs_DI_rm.eigenvalue",header=FALSE)
eigenvectors <- read.table("44_dataset_hapmap_autosomalSNPs_DI_rm.eigenvector")
# outliers <- paste( read.table("28_ceu_outliers",header=F)$V1,
#                    read.table("28_ceu_outliers",header=F)$V2,
#                    sep=":")

# DO NOT ENTER THE MERGED DATASET; ONLY DATASET WITH NON-REFERENCE INDIVIDUALS
sample_individuals <- paste( read.table("35_dataset_flipped_snps.fam")$V1,
                             read.table("35_dataset_flipped_snps.fam")$V2,
                             sep=":")

hapmap_individuals <- paste(read.table("/Users/naim/Documents/Strug/hapmap_data/07_hapmap_hg19_related_rm_ceu_tsi.fam")$V1,
                            read.table("/Users/naim/Documents/Strug/hapmap_data/07_hapmap_hg19_related_rm_ceu_tsi.fam")$V2,
                            sep=":")
relationships <- read.table("hapmap3r3_b36_relationship.txt", header=TRUE)
relationships$ID <- paste(relationships$FID, relationships$IID, sep=":")

#### SCREEPLOT
#png("20_screeplot.png")
plot(eigenvalues$V1, xlim=c(1,10), ylab="Eigenvalue",xlab="Component number", main="Screeplot") #screeplot
#dev.off()


#### Draw box
library(rgl)

# subset ethnic backgrounds and sample individuals
eigen_sample <- subset(eigenvectors, eigenvectors$V1 %in% sample_individuals)
eigen_hapmap <- subset(eigenvectors, eigenvectors$V1 %in% relationships$ID)
#eigen_outliers <- subset(eigenvectors, eigenvectors$V1 %in% outliers)
#eigen_outliers_5iter <- subset(eigenvectors, eigenvectors$V1 %in% outliers[1:15])
#eigen_outliers_20iter <- subset(eigenvectors, eigenvectors$V1 %in% outliers[16:32])
populations <- levels(relationships$population)
colours <- c("red", "green", "purple", "black", "brown", "blue", 
             "cadetblue", "pink", 54, 53, 55, "orange")

#draw the sample and hapmap individuals' first 3 PC points
with(eigenvectors, plot3d(V2,V3,V4, type="n"))
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_hapmap, eigen_hapmap$V1 %in% pop_i$ID)
  with(eigen_i, points3d(V2,V3,V4, col=colours[i+1]))
}
with(eigen_sample, points3d(V2,V3,V4,col=colours[1], size=6))
# eigen_sam_K <- subset(eigen_sample, grepl("(^[C,R][K,k].*)",as.character(eigen_sample$V1)))
# eigen_sam_S <- subset(eigen_sample, grepl("(^S.*)",as.character(eigen_sample$V1)))
# eigen_sam_rest <- subset(eigen_sample, !(as.character(eigen_sample$V1) %in% c(as.character(eigen_sam_K$V1),as.character(eigen_sam_S$V1))))
# with(eigen_sam_K, points3d(V2,V3,V4,col=colours[1],size=2))
# with(eigen_sam_S, points3d(V2,V3,V4,col="firebrick4",size=4))
# with(eigen_sam_rest, points3d(V2,V3,V4,col="skyblue",size=4))

#with(eigen_outliers, points3d(V2,V3,V4,col="firebrick4", size=10, pch="*"))
# with(eigen_outliers_5iter, points3d(V2,V3,V4,col="firebrick4", size=10, pch="*"))
# with(eigen_outliers_20iter, points3d(V2,V3,V4,col="orange", size=10, pch="*"))
# unknown_eth_ids <- paste(subset(ceu_tsi_sample_ind, ceu_tsi_sample_ind$Ethnicity %in% "Unknown")$FID,
#                          subset(ceu_tsi_sample_ind, ceu_tsi_sample_ind$Ethnicity %in% "Unknown")$IID,
#                          sep=":")
# eigen_unknown_eth <- subset(eigenvectors, eigenvectors$V1 %in% unknown_eth_ids)
#with(eigen_unknown_eth, points3d(V2,V3,V4,col="firebrick4", size=10, pch="*"))

# eigen_5iter_hapmap_coord <- subset(eigenvectors,eigenvectors$V1 %in% outliers[1:15])
# eigen_20iter_hapmap_coord <- subset(eigenvectors,eigenvectors$V1 %in% outliers[16:32])
# with(eigen_5iter_hapmap_coord, points3d(V2,V3,V4,col="firebrick4", size=10, pch="*"))
# with(eigen_20iter_hapmap_coord, points3d(V2,V3,V4,col="orange", size=10, pch="*"))


# Read in outliers:
outliers <- paste(read.table("46_outliers.txt",header=F)$V1,
                  read.table("46_outliers.txt",header=F)$V2,
                  sep=":")
outliers_removed <- subset(sample_individuals, !(sample_individuals %in% outliers))
#==========================================================================================================================
vectorIIDsTOlist <- function(samples) {
  fid <- iid <- NULL
  for (i in 1:length(samples)) {
    fid <- c(fid, gsub("(.*):.*","\\1",samples[i]))
    iid <- c(iid, gsub(".*:(.*)","\\1",samples[i]))
  }
  return(data.frame(FID=fid,IID=iid))
}
#==========================================================================================================================
# outliers_removed <- vectorIIDsTOlist(outliers_removed)
# write.table(outliers_removed, "relateds23_9CEU_list.txt", quote=F, col.names=F, row.names=F)
# 
# all_independent_ceu_cts_cases<-read.table("Step30_ceu_CTS_cases_only.fam",header=F)[,1:2]
# names(all_independent_ceu_cts_cases)<-c("FID","IID")
# final_list <- rbind(all_independent_ceu_cts_cases,outliers_removed)
# relationship_summary2(as.character(outliers_removed[,2]), rbind(parent_offsprings,sib_pairs,second_degree_pairs))
# temp<-relationship_summary2(as.character(all_independent_ceu_cts_cases[,2]), rbind(parent_offsprings,sib_pairs,second_degree_pairs))
# temp[temp$Related==TRUE,]
# temp<-relationship_summary2(as.character(final_list[,2]), rbind(parent_offsprings,sib_pairs,second_degree_pairs))
# final_list_relationships<-temp[temp$Related==TRUE,]
# setwd("/Users/naim/Desktop/")
# write.csv(outliers_removed,"temp.csv",row.names=F,quote=F)

# Draw the legend
par(mfrow=c(1,1))
plot(c(0,1),c(0,12),type='n',xlab='',ylab='')
x1<-rep(0,12)
x2<-rep(0.2,12)
y1=0:11
y2=1:12
rect(x1,y1,x2,y2,col=colours)
#x.text=rep(0.4,12)
#y.text=1:12
#text(x.text,y.text,express
text(0.4,0.5,expression('Samples'))
text(0.4,1.5,expression('ASW-African'))
text(0.4,2.5,expression('CEU-Caucasian'))
text(0.4,3.5,expression('CHB-Chinese-Beijing'))
text(0.4,4.5,expression('CHD-Chinese-Denver'))
text(0.4,5.5,expression('GIH-Indian'))
text(0.4,6.5,expression('JPT-Japanese'))
text(0.4,7.5,expression('LWK-Luhya'))
text(0.4,8.5,expression('MEX-Mexican'))
text(0.4,9.5,expression('MKK-Maasia'))
text(0.4,10.5,expression('TSI-Italian'))
text(0.4,11.5,expression('YRI-Yoruba'))


# Want to determine ethnicity of each individual
pop_i_samples <- list()
esample_individuals <- eigen_sample
print(paste("Number of individuals per population"))
populations <- levels(relationships$population)
for (i in 1:length(populations)){
  pop_i <- subset(relationships, relationships$population %in% populations[i])
  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(eigen_hapmap, eigen_hapmap$V1 %in% pop_i$ID)
  individual_indices <- isInBox(eigen_i[,2:4],esample_individuals[,2:4]) 
  print(paste(populations[i],sum(individual_indices)))
  pop_i_samples[[i]] <- as.character(esample_individuals$V1[individual_indices])
  names(pop_i_samples)[i] <- populations[i]
  if (sum(individual_indices)>0) esample_individuals <- esample_individuals[-which(individual_indices),] # remove sample individuals to prevent duplicates
}
sapply(pop_i_samples,length)
sum(sapply(pop_i_samples,length))

# who has not been classified?
esample_individuals$V1
dim(esample_individuals)
with(esample_individuals, points3d(V2,V3,V4, col="grey", size=8, pch="*"))

# include last category of unclassified individuals:
pop_i_samples[[12]] <- as.character(esample_individuals$V1) 
populations <- c(populations, "Unknown")
names(pop_i_samples)[12] <- populations[12]


# Convert pop_i_samples to a data.frame:
ethnicity<-NULL
fid<-NULL
iid<-NULL
for (i in 1:length(pop_i_samples)) {
  ethnicity <- c(ethnicity,rep(names(pop_i_samples)[i], length(pop_i_samples[[i]])))
  fid <- c(fid, gsub("(.*):.*","\\1",pop_i_samples[[i]]))
  iid <- c(iid, gsub(".*:(.*)","\\1",pop_i_samples[[i]]))
}
ind_ethnicities <- data.frame(FID=fid, IID=iid, Ethnicity=ethnicity)
write.table(ind_ethnicities, "ethnicities.txt",
            quote=F, row.names=F, col.names=T)



getFID <- function(x) return(gsub("(.*):.*","\\1", x))
getIID <- function(x) return(gsub(".*:(.*)","\\1", x))



### List the sample individuals plus Hapmap individuals to run
### smartpca on to check for outliers more than 6 stdevs away

# Outliers away from the CEU/TSI cluster
# 1. Remove sample individuals that fall into other clusters:
ceu_tsi_sample_ind <- subset(ind_ethnicities, ind_ethnicities$Ethnicity %in% c("CEU", "TSI", "Unknown"))
dim(ceu_tsi_sample_ind)
# 134

# 2. List Hapmap CEU and TSI individuals only:
ceu_tsi_hapmap_ind <- subset(relationships, relationships$population %in% c("CEU","TSI"))[,c(8,7)]
hapmap_ceu_tsi_fids <- as.character(sapply(ceu_tsi_hapmap_ind$ID, getFID))
hapmap_ceu_tsi_iids <- as.character(sapply(ceu_tsi_hapmap_ind$ID, getIID))
ceu_tsi_hapmap_ind <- cbind(hapmap_ceu_tsi_fids, hapmap_ceu_tsi_iids, as.character(ceu_tsi_hapmap_ind[,2]))
colnames(ceu_tsi_hapmap_ind) <- c("FID","IID","Ethnicity")
dim(ceu_tsi_hapmap_ind)
# 267

# 3. Merge them for extraction using PLINK and subsequent PCA to detect individuals 6 stdevs away
ceu_merge <- rbind(ceu_tsi_sample_ind, ceu_tsi_hapmap_ind)
ceu_ind_for_pca <- write.table(ceu_merge, "ceu_ind_for_pca.txt",
                               quote=F, row.names=F, col.names=F)




# 4. After identifying the outliers, create final CEU/TSI list of individuals to keep and analyze
ceu_tsi_sample_ids <- paste(ceu_tsi_sample_ind$FID,ceu_tsi_sample_ind$IID,sep=":")
final_ceu <- subset(ceu_tsi_sample_ind, !(ceu_tsi_sample_ids %in% outliers[1:15]))
unknown_CEU_ind <- subset(final_ceu,final_ceu$Ethnicity %in% "Unknown")
# 26 Unknown individuals classified as CEU

final_ceu$Ethnicity <- as.character(final_ceu$Ethnicity)
final_ceu$Ethnicity[which(final_ceu$Ethnicity %in% "Unknown")] <- "CEU-TSI"

write.table(final_ceu, "ceu_individuals.txt",
            quote=F, row.names=F, col.names=F)

###########################
# What about the related pairs that were removed for PCA?
###########################
# Types of relatedness: siblings, parent-child, 2nd-degree, twins

pop_relatives <- list()
siblings_list <- c(as.character(sib_pairs$IID1), 
                   as.character(sib_pairs$IID2))
parents_list <- c(as.character(parent_offspring$IID1), 
                  as.character(parent_offspring$IID2))
twins_list <- c(as.character(twins_duplicates$IID1),
                as.character(twins_duplicates$IID2))
second_degree_list <- c(as.character(second_degree_relatives$IID1), 
                        as.character(second_degree_relatives$IID2))
matching_sibs = matching_parents = matching_twins = matching_second_degree = NULL
for (i in 1:length(pop_i_samples)) { # for each population
  pop_i <- data.frame(ID=pop_i_samples[[i]])
  if (dim(pop_i)[1] > 0) { # if there are any individuals in this population
    pop_i$FID <- sapply(strsplit(as.character(pop_i$ID),":"),"[[",1)
    pop_i$IID <- sapply(strsplit(as.character(pop_i$ID),":"),"[[",2)
    pop_i <- pop_i[,2:3]
    pop_i$siblings <- 0
    pop_i$parents <- 0
    pop_i$twins <- 0
    pop_i$second_degree <- 0
    
    for (j in 1:dim(pop_i)[1]) { # for each individual in population i
      # get how many siblings
      pop_i$siblings[j] <- sum(siblings_list %in% as.character(pop_i$IID[j]))
      if (length(which(as.character(sib_pairs$IID1) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_sibs <- c(matching_sibs, 
                           as.character(sib_pairs$IID2[which(as.character(sib_pairs$IID1) 
                                                             %in% as.character(pop_i$IID[j]))]))
      } 
      if (length(which(as.character(sib_pairs$IID2) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_sibs <- c(matching_sibs, 
                           as.character(sib_pairs$IID1[which(as.character(sib_pairs$IID2) 
                                                             %in% as.character(pop_i$IID[j]))]))
      }
      
      # how many parents
      pop_i$parents[j] <- sum(parents_list %in% as.character(pop_i$IID[j]))
      if (length(which(as.character(parent_offspring$IID1) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_parents <- c(matching_parents, 
                              as.character(parent_offspring$IID2[which(as.character(parent_offspring$IID1) 
                                                                       %in% as.character(pop_i$IID[j]))]))
      } 
      if (length(which(as.character(parent_offspring$IID2) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_parents <- c(matching_parents, 
                              as.character(parent_offspring$IID1[which(as.character(parent_offspring$IID2) 
                                                                       %in% as.character(pop_i$IID[j]))]))
      }
      
      
      # how many twins/duplicates
      pop_i$twins[j] <- sum(twins_list %in% as.character(pop_i$IID[j]))
      if (length(which(as.character(twins_duplicates$IID1) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_twins <- c(matching_twins, 
                            as.character(twins_duplicates$IID2[which(as.character(twins_duplicates$IID1) 
                                                                     %in% as.character(pop_i$IID[j]))]))
      } 
      if (length(which(as.character(twins_duplicates$IID2) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_twins <- c(matching_twins, 
                            as.character(twins_duplicates$IID1[which(as.character(twins_duplicates$IID2) 
                                                                     %in% as.character(pop_i$IID[j]))]))
      }
      
      # how many 2nd degree relatives
      pop_i$second_degree[j] <- sum(second_degree_list%in% as.character(pop_i$IID[j]))
      if (length(which(as.character(second_degree_relatives$IID1) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_second_degree <- c(matching_second_degree, 
                                    as.character(second_degree_relatives$IID2[which(as.character(second_degree_relatives$IID1) 
                                                                                    %in% as.character(pop_i$IID[j]))]))
      } 
      if (length(which(as.character(second_degree_relatives$IID2) %in% as.character(pop_i$IID[j]))) > 0) 
      {
        matching_second_degree <- c(matching_second_degree, 
                                    as.character(second_degree_relatives$IID1[which(as.character(second_degree_relatives$IID2) 
                                                                                    %in% as.character(pop_i$IID[j]))]))
      }
      
    }
  }
  pop_relatives[[i]] <- pop_i
  names(pop_relatives)[i] <- populations[i]
}

pop_relatives

summary_table <- NULL
j<-0
for (i in 1:length(pop_relatives)) {
  if (dim(pop_relatives[[i]])[1]>0) {
    summary_table <- rbind(summary_table, sapply(pop_relatives[[i]][,3:6],sum))
    j<-j+1
    row.names(summary_table)[j] <- populations[i]
  }
}

summary_table

# number of independent individuals in each population
(indep_summary_table <- sapply(pop_i_samples,length))
# number of relatives per population
(relatives_summary_tables <- apply(summary_table,1,sum))

(new_indep_summary_table <- indep_summary_table[names(relatives_summary_tables)])

(totals <- new_indep_summary_table + relatives_summary_tables)
sum(totals)


#==============================
all_relatives_ids <- c(matching_parents,matching_second_degree,matching_sibs,matching_twins)
length(all_relatives_ids)
#461
length(unique(all_relatives_ids))
#395


#============================================
# Remove MKK and MEX populations from Hapmap and redo PCA 
#============================================
new_populations <- populations[-which(populations %in% c("MKK","MEX"))]
new_relationships <- subset(relationships, relationships$population %in% new_populations)

new_eigen_hapmap <- subset(eigenvectors, eigenvectors$V1 %in% new_relationships$ID)
new_eigenvectors <- rbind(new_eigen_hapmap, eigen_sample)
new_colours <- c("red", "green", "black", "brown", "blue", 
                 "cadetblue", "pink", 54, 55, "orange")

ind_to_remove <- subset(relationships, relationships$population %in% c("MKK","MEX"))
ind_to_remove <- ind_to_remove[,1:2]
write.table(ind_to_remove, "MKK_MEX_hapmap_individuals_to_remove.txt", 
            quote=FALSE,row.names=FALSE, col.names=FALSE)

# Draw box
library(rgl)
with(new_eigenvectors, plot3d(V2,V3,V4, type="n"))

#draw the sample and hapmap individuals' first 3 PC points
with(eigen_sample, points3d(V2,V3,V4,col=new_colours[1], size=6))
for (i in 1:length(new_populations)){
  pop_i <- subset(new_relationships, new_relationships$population %in% new_populations[i])
  pop_i$ID <- paste(pop_i$FID,pop_i$IID,sep=":")
  eigen_i <- subset(new_eigen_hapmap, new_eigen_hapmap$V1 %in% pop_i$ID)
  with(eigen_i, points3d(V2,V3,V4, col=new_colours[i+1]))
}

# Draw the legend
par(mfrow=c(1,1))
plot(c(0,1),c(0,10),type='n',xlab='',ylab='')
x1<-rep(0,9)
x2<-rep(0.2,10)
y1=0:9
y2=1:10
rect(x1,y1,x2,y2,col=new_colours)
#x.text=rep(0.4,12)
#y.text=1:12
#text(x.text,y.text,express
text(0.4,0.5,expression('UKRE Data'))
text(0.4,1.5,expression('ASW-African'))
text(0.4,2.5,expression('CEU-Caucasian'))
text(0.4,3.5,expression('CHB-Chinese-Beijing'))
text(0.4,4.5,expression('CHD-Chinese-Denver'))
text(0.4,5.5,expression('GIH-Indian'))
text(0.4,6.5,expression('JPT-Japanese'))
text(0.4,7.5,expression('LWK-Luhya'))
text(0.4,8.5,expression('TSI-Italian'))
text(0.4,9.5,expression('YRI-Yoruba'))


