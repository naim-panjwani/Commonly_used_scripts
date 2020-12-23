# R functions useful for association analyses

########## Global variables ###########
unaffected.code<-1
affected.code<-2
male.sex.code<-1
female.sex.code<-2
########### FUNCTIONS #################
removeZeros <- function(allele1, allele2) {
  newAllele1<-allele1
  newAllele2<-allele2
  if (any(newAllele1 %in% 0)) {
    pos<-which(newAllele1 %in% 0)
    newAllele1<-factor(newAllele1[-(pos)])
    newAllele2<-factor(newAllele2[-(pos)]) }
  if (any(newAllele2 %in% 0)) {
    pos<-which(newAllele2 %in% 0)
    newAllele1<-factor(newAllele1[-(pos)])
    newAllele2<-factor(newAllele2[-(pos)]) } 
  return(list(A1=newAllele1,A2=newAllele2))
}
MinorAllele <- function(allele1, allele2) {
  # Given vectors of allele1 and allele2 (each row is per individual),
  # return the allele nucleotide letter in lowest proportion 
  if (length(allele1) != length(allele2)) {
    return (c("Allele vector lengths differ"))
    } else {
      newAlleleList<-removeZeros(allele1,allele2)
      allele1<-newAlleleList$A1
      allele2<-newAlleleList$A2 
      tempF<-factor(cbind(as.character(allele1),as.character(allele2))) 
      if (length(levels(tempF)) < 2) {
        c("Only one allele type")
      } else {
        y<-character(1)
        if (table(tempF)[1] <= table(tempF)[2]) {
          y <- levels(tempF)[1]
          } else {
            y <- levels(tempF)[2] }
    return(y) }
    } 
}
MajorAllele <- function(allele1, allele2) {
  # Given vectors of allele1 and allele2 (each row is per individual),
  # return the allele nucleotide letter in highest proportion 
  if (length(allele1) != length(allele2)) {
    return (c("Allele vector lengths differ"))
    } else {
      newAlleleList<-removeZeros(allele1,allele2)
      allele1<-newAlleleList$A1
      allele2<-newAlleleList$A2
      tempF<-factor(cbind(as.character(allele1),as.character(allele2))) 
      if (length(levels(tempF)) < 2) {
        c("Only one allele type")
        } else {
          y<-character(1)
          if (table(tempF)[1] >= table(tempF)[2]) {
            y <- levels(tempF)[1]
            } else {
              y <- levels(tempF)[2] }
          return(y) }
    }
}
GenCounts <- function(allele1, allele2) {
  # Given allele1, allele2 vectors of a SNP
  # returns the genotype counts for each AA, Aa, aa types 
  newAlleleList<-removeZeros(allele1,allele2)
  allele1<-newAlleleList$A1
  allele2<-newAlleleList$A2
  minor <- MinorAllele(allele1, allele2)
  major <- MajorAllele(allele1, allele2)
  AA.count<-0
  Aa.count<-0
  aa.count<-0
  unknown.count<-0
  error<-matrix(ncol=2)
  if (minor %in% major) {
    c("SNP is not heterozygous at this site") } else {
      for (i in 1:length(allele1)) { #for each individual
        if ((allele1[i] == major) & allele2[i] == major) { #if aa
          aa.count<-aa.count+1
        } else if ((allele1[i] == major & allele2[i] == minor) |
                     (allele1[i] == minor & allele2[i] == major)) { #if Aa or aA
          Aa.count<-Aa.count+1
        } else if (allele1[i] == minor & allele2[i] == minor) { #if AA AA.count<-AA.count+1
        } else { #for troubleshooting purposes unknown.count<-unknown.count+1 error<-rbind(error,c(allele1[i],allele2[i]))
        } }
      return(list(result=c(AA.count,Aa.count,aa.count,unknown.count), error=error[-1,]))
    } }
GenCountsMajor <- function(allele1, allele2) {
  # Given allele1, allele2 vectors of a SNP
  # returns the genotype counts for each AA, Aa, aa types 
  newAlleleList<-removeZeros(allele1,allele2) 
  allele1<-newAlleleList$A1
  allele2<-newAlleleList$A2
  minor <- MinorAllele(allele1, allele2)
  major <- MajorAllele(allele1, allele2)
  AA.count<-0
  Aa.count<-0
  aa.count<-0
  unknown.count<-0
  error<-matrix(ncol=2)
  if (minor %in% major) {
    c("SNP is not heterozygous at this site")
    } else {
      for (i in 1:length(allele1)) { #for each individual
        if ((allele1[i] == major) & allele2[i] == major) { #if aa
          AA.count<-AA.count+1
        } else if ((allele1[i] == major & allele2[i] == minor) |
                     (allele1[i] == minor & allele2[i] == major)) { #if Aa or aA
          Aa.count<-Aa.count+1
        } else if (allele1[i] == minor & allele2[i] == minor) { #if AA
          aa.count<-aa.count+1
        } else { #for troubleshooting purposes
          unknown.count<-unknown.count+1
          error<-rbind(error,c(allele1[i],allele2[i])) 
        }
      }
      return(list(result=c(AA.count,Aa.count,aa.count,unknown.count), error=error[-1,]))
    }
}


Genotype <- function(allele1, allele2, AA_code=2, Aa_code=1, aa_code=0) {
  # Given allele1 and allele2 vectors of a SNP,
  # returns genotype type "AA", "Aa" or "aa" and
  # allele coding in data.frame format
  # MISSING ALLELES ARE ALLOWED AND WON'T BE REMOVED
  newAllele1<-allele1
  newAllele2<-allele2
  minor<-MinorAllele(newAllele1,newAllele2)
  major<-MajorAllele(newAllele1,newAllele2)
  result<-data.frame(ncol=2,nrow=length(newAllele1))
  colnames(result)<-c("Genotype", "Genotype.code")
  for(i in 1:length(newAllele1)) {
    wildtype<-newAllele1[i] == major & newAllele2[i] == major
    het<-(newAllele1[i] == major & newAllele2[i] == minor) |
      (newAllele1[i] == minor & newAllele2[i] == major)
    hom<-newAllele1[i] == minor & newAllele2[i] == minor
    if (wildtype) { #if aa
        result[i,] <- c("aa", as.numeric(aa_code))
      } else if (het) { #if Aa or aA
        result[i,] <- c("Aa", as.numeric(Aa_code))
      } else if (hom) { #if AA
        result[i,] <- c("AA", as.numeric(AA_code))
      } else { #for troubleshooting purposes
        result[i,] <- c("Missing", NA)
      }
  }
  return(result)
}

GenotypeMajor <- function(allele1, allele2, AA_code=2, Aa_code=1, aa_code=0) {
  # Given allele1 and allele2 vectors of a SNP,
  # returns genotype type "AA", "Aa" or "aa" and
  # allele coding in data.frame format
  # MISSING ALLELES ARE ALLOWED AND WON'T BE REMOVED 
  newAllele1<-allele1
  newAllele2<-allele2
  minor<-MinorAllele(newAllele1,newAllele2) 
  major<-MajorAllele(newAllele1,newAllele2) 
  result<-data.frame(ncol=2,nrow=length(newAllele1)) 
  colnames(result)<-c("Genotype", "Genotype.code") 
  for(i in 1:length(newAllele1)) {
    wildtype<-newAllele1[i] == minor & newAllele2[i] == minor 
    het<-(newAllele1[i] == major & newAllele2[i] == minor) |
                 (newAllele1[i] == minor & newAllele2[i] == major) 
    hom<-newAllele1[i] == major & newAllele2[i] == major 
    if (wildtype) { #if aa
      result[i,] <- c("aa", as.numeric(aa_code)) 
      } else if (het) { #if Aa or aA
        result[i,] <- c("Aa", as.numeric(Aa_code))
        } else if (hom) { #if AA
          result[i,] <- c("AA", as.numeric(AA_code))
          } else { #for troubleshooting purposes
            result[i,] <- c("Missing", NA)
          }
  }
  return(result)
}


SNPDatum <- function(ped, allele1, allele2,
                     AA_code=2, Aa_code=1, aa_code=0) { #default; change for DOMDEV
  # Given the first 6 columns of the pedigree data file,
  # returns a dataframe with
  # FID, IID, Sex code, Phenotype, Genotype, Genotype Code
  result<-data.frame(FID=ped[,1], IID=ped[,2], Sex=ped[,5], Phenotype=ped[,6], 
                     Geno=Genotype(allele1,allele2, AA_code, Aa_code, aa_code)[,1],
                     Geno.code=Genotype(allele1,allele2, as.numeric(AA_code), as.numeric(Aa_code), as.numeric(aa_code))[,2]) 
  minor <- MinorAllele(allele1, allele2)
  major <- MajorAllele(allele1, allele2) 
  newAllele1<-allele1
  newAllele2<-allele2
  #removing missing data
  if (any(is.na(result[,"Geno.code"]))) {
    pos<-which(is.na(result[,"Geno.code"]))
    result<-result[-(pos),] }
  if (any(is.na(result[,"Geno.code"]))) {
    pos<-which(is.na(result[,"Geno.code"]))
    result<-result[-(pos),]
  }
  return(result)
}

SNPDatumMajor <- function(ped, allele1, allele2,
                          AA_code=2, Aa_code=1, aa_code=0) { #default; change for DOMDEV
  # Given the first 6 columns of the pedigree data file,
  # returns a dataframe with
  # FID, IID, Sex code, Phenotype, Genotype, Genotype Code 
  result<-data.frame(FID=ped[,1], IID=ped[,2], Sex=ped[,5], Phenotype=ped[,6], 
                     Geno=GenotypeMajor(allele1,allele2, AA_code, Aa_code, aa_code)[,1], 
                     Geno.code=GenotypeMajor(allele1,allele2, as.numeric(AA_code), as.numeric(Aa_code), as.numeric(aa_code))[,2]) 
  minor <- MinorAllele(allele1, allele2)
  major <- MajorAllele(allele1, allele2) 
  newAllele1<-allele1
  newAllele2<-allele2
  #removing missing data
  if (any(is.na(result[,"Geno.code"]))) {
    pos<-which(is.na(result[,"Geno.code"]))
    result<-result[-(pos),] }
  if (any(is.na(result[,"Geno.code"]))) {
    pos<-which(is.na(result[,"Geno.code"]))
    result<-result[-(pos),]
  }
  return(result)
}


Table23 <- function(pedfile, mapfile) {
  # Given a pedfile and mapfile, returns a 2x3 genotypic table of: 
  # $Genotypes: all case and control genotypes
  # $Controls: controls only
  # $Cases: cases only
  num.snps<-dim(mapfile)[1]
  pedCtrls <- subset(pedfile, V6 %in% 1)
  numCtrls <- dim(pedCtrls)[1]
  pedCases <- subset(pedfile, V6 %in% 2)
  numCases <- dim(pedCases)[1]
  ctrl.gen.counts<-matrix(nrow=num.snps,ncol=4) 
  colnames(ctrl.gen.counts) <- c("AA","Aa","aa","Errors") 
  rownames(ctrl.gen.counts) <- paste("CHR", mapfile$V1, mapfile$V2, sep=" ")
  case.gen.counts<-matrix(nrow=num.snps,ncol=4) 
  colnames(case.gen.counts) <- c("AA","Aa","aa","Errors") 
  rownames(case.gen.counts) <- rownames(ctrl.gen.counts) <- paste("CHR", mapfile$V1, mapfile$V2, sep=" ")
  
  for (i in 1:num.snps) { #for each SNP
    col.position <- 7 + 2*(i-1)
    for (j in 1:numCtrls) { # for each control individual
      ctrl.gen.counts[i,] <- GenCounts(pedCtrls[,col.position],pedCtrls[,col.position+1])$result
    }
    for (j in 1:numCases) { # for each case individual
      case.gen.counts[i,] <- GenCounts(pedCases[,col.position],pedCases[,col.position+1])$result
    }
  }
  control.genotypes <- ctrl.gen.counts[,-4] 
  case.genotypes <- case.gen.counts[,-4] 
  all.genotypes<-matrix(nrow=2*num.snps, ncol=3) 
  rnames<-character(2*num.snps)
  for (i in 1:num.snps) {
    all.genotypes[2*i-1,] <- control.genotypes[i,] 
    rnames[2*i-1] <- paste("CHR",mapfile$V1[i],mapfile$V2[i],"Controls",sep=".") 
    all.genotypes[2*i,] <- case.genotypes[i,]
    rnames[2*i] <- paste("CHR",mapfile$V1[i],mapfile$V2[i],"Cases",sep=".") 
  }
  rownames(all.genotypes) <- rnames
  colnames(all.genotypes) <- c("AA","Aa","aa")
  return(list(Genotypes=all.genotypes, Controls=control.genotypes, Cases=case.genotypes))
}


Table23Major <- function(pedfile, mapfile) {
  # Given a pedfile and mapfile, returns a 2x3 genotypic table of: 
  # $Genotypes: all case and control genotypes
  # $Controls: controls only
  # $Cases: cases only
  num.snps<-dim(mapfile)[1]
  pedCtrls <- subset(pedfile, V6 %in% 1)
  numCtrls <- dim(pedCtrls)[1]
  pedCases <- subset(pedfile, V6 %in% 2)
  numCases <- dim(pedCases)[1]
  ctrl.gen.counts<-matrix(nrow=num.snps,ncol=4) 
  colnames(ctrl.gen.counts) <- c("AA","Aa","aa","Errors") 
  rownames(ctrl.gen.counts) <- paste("CHR", mapfile$V1, mapfile$V2, sep=" ")
  case.gen.counts<-matrix(nrow=num.snps,ncol=4) 
  colnames(case.gen.counts) <- c("AA","Aa","aa","Errors") 
  rownames(case.gen.counts) <- rownames(ctrl.gen.counts) <-
    paste("CHR", mapfile$V1, mapfile$V2, sep=" ")
  for (i in 1:num.snps) { #for each SNP
    col.position <- 7 + 2*(i-1)
    for (j in 1:numCtrls) { # for each control individual
      ctrl.gen.counts[i,] <- GenCountsMajor(pedCtrls[,col.position],pedCtrls[,col.position+1])$result
    }
    for (j in 1:numCases) { # for each case individual
      case.gen.counts[i,] <- GenCountsMajor(pedCases[,col.position],pedCases[,col.position+1])$result
    } }
  control.genotypes <- ctrl.gen.counts[,-4] 
  case.genotypes <- case.gen.counts[,-4] 
  all.genotypes<-matrix(nrow=2*num.snps, ncol=3) 
  rnames<-character(2*num.snps)
  for (i in 1:num.snps) {
    all.genotypes[2*i-1,] <- control.genotypes[i,] 
    rnames[2*i-1]<-
      paste("CHR",mapfile$V1[i],mapfile$V2[i],"Controls",sep=".") 
    all.genotypes[2*i,] <- case.genotypes[i,]
    rnames[2*i] <-
      paste("CHR",mapfile$V1[i],mapfile$V2[i],"Cases",sep=".") }
  rownames(all.genotypes) <- rnames
  colnames(all.genotypes) <- c("AA","Aa","aa")
  return(list(Genotypes=all.genotypes, Controls=control.genotypes, Cases=case.genotypes))
}

getMAF <- function(allele1, allele2) {
  temp <- as.numeric(Genotype(allele1,allele2, AA_code=2, Aa_code=1, aa_code=0)$Genotype.code)
  tot <- 2*(length(temp))
  A.count <- sum(temp)
  maf <- A.count / tot
  return(maf)
}

######################################################################
############# General Logistic Function ############
# Input: pedfile, mapfile, model: multiplicative/additive (default), genotypic, recessive, dominant
# Output: $REG: a matrix of 1 SNP per row; columns: b1 OR, b1 pvalue, b2 OR, b2 pvalue, AIC
# $Table: the corresponding genotypic table
######################################################################

Logistic <- function(pedfile, mapfile, model = "multiplicative", sex=TRUE) {
  num.snps<-dim(mapfile)[1]
  pedCtrls <- subset(pedfile, V6 %in% 1)
  numCtrls <- dim(pedCtrls)[1]
  pedCases <- subset(pedfile, V6 %in% 2)
  numCases <- dim(pedCases)[1]
  AA_code<-2; Aa_code<-1; aa_code<-0
  AA_codeD<-0; Aa_codeD<-1; aa_codeD<-0
  logistic<-list()
  length(logistic)<-num.snps
  snp.data=list()
  length(snp.data)<-num.snps
  domdev.data<-list()
  length(domdev.data)<-num.snps
  if (model %in% "multiplicative" | model %in% "additive") { AA_code=2; Aa_code=1; aa_code=0
  } else if (model == "genotypic") {
    AA_code=2; Aa_code=1; aa_code=0
    AA_codeD<-0; Aa_code<-1; aa_code<-0
  } else if (model == "recessive") {
    AA_code=1; Aa_code=0; aa_code=0
  } else if (model == "dominant") {
    AA_code=1; Aa_code=1; aa_code=0
  }
  Logistic.results<-matrix()
  b1.pvalues<-numeric(num.snps)
  b2.pvalues<-numeric(num.snps)
  for (i in 1:num.snps) { #for each SNP
    col.position <- 7 + 2*(i-1)
    snp.data[[i]]<-SNPDatum(pedfile[,1:6],
                            pedfile[,col.position],
                            pedfile[,col.position+1],
                            AA_code, Aa_code, aa_code)
    domdev.data[[i]]<-SNPDatum(pedfile[,1:6], pedfile[,col.position],
                               pedfile[,col.position+1],
                               AA_codeD, Aa_codeD, aa_codeD)
    y <- snp.data[[i]]$Phenotype - 1 #Unaffected/controls are
    reference
    b1.ADD <- as.numeric(as.character(snp.data[[i]]$Geno.code))
    #wildtype (aa) is reference
    sex.covar <- snp.data[[i]]$Sex - 1 #males are reference
    b2.DOMDEV <- as.numeric(as.character(domdev.data[[i]]$Geno.code))
    if (model == "genotypic") {
      if (sex) {
        logistic[[i]]<-list(REG=glm(y~b1.ADD+b2.DOMDEV+sex.covar, family=binomial))
      } else { logistic[[i]]<-list(REG=glm(y~b1.ADD+b2.DOMDEV,
                                           family=binomial))
      }
    } else {
      if (sex) {
        logistic[[i]]<-list(REG=glm(y~b1.ADD+sex.covar, family=binomial))
      } else {
        logistic[[i]]<-list(REG=glm(y~b1.ADD, family=binomial))
      } }
    if (i==1) { Logistic.results<-t(sapply(logistic[[i]],function(x) exp(coef(x)))) 
                Logistic.results<-cbind(Logistic.results, AIC=logistic[[i]]$REG$aic)
    } else {
      nextrow<-c(t(sapply(logistic[[i]],function(x) exp(coef(x)))), AIC=logistic[[i]]$REG$aic)
      Logistic.results<-rbind(Logistic.results, nextrow) }
    if (model == "genotypic") {
      if (!is.na(Logistic.results[i,"b1.ADD"])){
        b1.pvalues[i]<- summary(logistic[[i]]$REG)$coefficients["b1.ADD","Pr(>|z|)"]
      } else {b1.pvalues[i]<-NA}
      if (!is.na(Logistic.results[i,"b2.DOMDEV"])) {
        b2.pvalues[i]<- summary(logistic[[i]]$REG)$coefficients["b2.DOMDEV","Pr(>|z|)"]
      } else {b2.pvalues[i]<-NA}
    } else {
      if (!is.na(Logistic.results[i,"b1.ADD"])){ 
        b1.pvalues[i] <- summary(logistic[[i]]$REG)$coefficients["b1.ADD","Pr(>|z|)"]
      } else {b1.pvalues[i]<-NA}
    } }
  if (model == "genotypic") { 
    final_logistic<-cbind(Logistic.results[,"b1.ADD"], b1.pvalues,
                          Logistic.results[,"b2.DOMDEV"], b2.pvalues,
                          Logistic.results[,"AIC"]) 
    colnames(final_logistic)<-c("b1.ADD.OR","b1.P", "b2.DOMDEV OR","b2.P", "AIC") 
    rownames(final_logistic)<-paste("CHR",mapfile$V1,mapfile$V2,sep=" ")
  } else {
    final_logistic<-cbind(Logistic.results[,"b1.ADD"], b1.pvalues, Logistic.results[,"AIC"])
    colnames(final_logistic)<-c("b1.ADD.OR","b1.P", "AIC")
    rownames(final_logistic)<-paste("CHR",mapfile$V1,mapfile$V2,sep=" ")
  }
  #### Return the correct 2x3 or 2x2 table ###### Gen.table<-matrix()
  if (model == "genotypic") {
    Gen.table<-Table23(pedfile,mapfile)
    } else if (model == "recessive") {
      tmp23<-Table23(pedfile,mapfile)$Genotypes 
      unaffecteds<-numeric(num.snps)
      for(i in 1:num.snps) {unaffecteds[i]<-tmp23[i,2]+tmp23[i,3]} 
      Gen.table<-cbind(tmp23[,1], unaffecteds) 
      colnames(Gen.table)<-c("AA", "Aa/aa")
    } else if (model == "dominant"){ 
      tmp23<-Table23(pedfile,mapfile)$Genotypes 
      affecteds<-numeric(num.snps)
      for(i in 1:num.snps) {affecteds[i]<-tmp23[i,1]+tmp23[i,2]} 
      Gen.table<-cbind(affecteds,tmp23[,3]) 
      colnames(Gen.table)<-c("AA/Aa","aa")
    } else { #additive/multiplicative Gen.table<-Table23(pedfile,mapfile)
    }
  return(list(REG=final_logistic, Table=Gen.table, LOG=logistic)) 
}


LogisticMajor <- function(pedfile, mapfile, model = "multiplicative", sex=TRUE) {
  num.snps<-dim(mapfile)[1]
  pedCtrls <- subset(pedfile, V6 %in% 1)
  numCtrls <- dim(pedCtrls)[1]
  pedCases <- subset(pedfile, V6 %in% 2)
  numCases <- dim(pedCases)[1]
  AA_code<-2; Aa_code<-1; aa_code<-0
  AA_codeD<-0; Aa_codeD<-1; aa_codeD<-0
  logistic<-list()
  length(logistic)<-num.snps
  snp.data=list()
  length(snp.data)<-num.snps
  domdev.data<-list()
  length(domdev.data)<-num.snps
  if (model %in% "multiplicative" | model %in% "additive") { AA_code=2; Aa_code=1; aa_code=0
  } else if (model == "genotypic") {
    AA_code=2; Aa_code=1; aa_code=0
    AA_codeD<-0; Aa_code<-1; aa_code<-0
  } else if (model == "recessive") {
    AA_code=1; Aa_code=0; aa_code=0
  } else if (model == "dominant") {
    AA_code=1; Aa_code=1; aa_code=0
  }
  Logistic.results<-matrix()
  b1.pvalues<-numeric(num.snps)
  b2.pvalues<-numeric(num.snps)
  for (i in 1:num.snps) { #for each SNP
    col.position <- 7 + 2*(i-1) 
    snp.data[[i]]<-SNPDatumMajor(pedfile[,1:6], pedfile[,col.position],pedfile[,col.position+1],AA_code, Aa_code, aa_code) 
    domdev.data[[i]]<-SNPDatumMajor(pedfile[,1:6],pedfile[,col.position],pedfile[,col.position+1],AA_codeD, Aa_codeD, aa_codeD)
    y <- snp.data[[i]]$Phenotype - 1 #Unaffected/controls are reference
    b1.ADD <- as.numeric(as.character(snp.data[[i]]$Geno.code)) #wildtype (aa) is reference
    sex.covar <- snp.data[[i]]$Sex - 1 #males are reference
    b2.DOMDEV <- as.numeric(as.character(domdev.data[[i]]$Geno.code))
    if (model == "genotypic") {
      if (sex) {
        logistic[[i]]<-list(REG=glm(y~b1.ADD+b2.DOMDEV+sex.covar, family=binomial))
      } else { logistic[[i]]<-list(REG=glm(y~b1.ADD+b2.DOMDEV,
                                           family=binomial))
      }
    } else {
      if (sex) {
        logistic[[i]]<-list(REG=glm(y~b1.ADD+sex.covar, family=binomial))
      } else {
        logistic[[i]]<-list(REG=glm(y~b1.ADD, family=binomial))
      } }
    if (i==1) { Logistic.results<-t(sapply(logistic[[i]],function(x) exp(coef(x)))) 
      Logistic.results<-cbind(Logistic.results, AIC=logistic[[i]]$REG$aic)
    } else {
      nextrow<-c(t(sapply(logistic[[i]],function(x) exp(coef(x)))), AIC=logistic[[i]]$REG$aic)
      Logistic.results<-rbind(Logistic.results, nextrow) }
    if (model == "genotypic") {
      if (!is.na(Logistic.results[i,"b1.ADD"])){
        b1.pvalues[i]<- summary(logistic[[i]]$REG)$coefficients["b1.ADD","Pr(>|z|)"]
      } else {b1.pvalues[i]<-NA}
      if (!is.na(Logistic.results[i,"b2.DOMDEV"])) {
        b2.pvalues[i]<- summary(logistic[[i]]$REG)$coefficients["b2.DOMDEV","Pr(>|z|)"]
      } else {b2.pvalues[i]<-NA}
    } else {
      if (!is.na(Logistic.results[i,"b1.ADD"])){
        b1.pvalues[i]<- summary(logistic[[i]]$REG)$coefficients["b1.ADD","Pr(>|z|)"]
      } else {b1.pvalues[i]<-NA}
    }
  }
  if (model == "genotypic") { 
    final_logistic<-cbind(Logistic.results[,"b1.ADD"], b1.pvalues,
                                                    Logistic.results[,"b2.DOMDEV"], b2.pvalues,
                                                    Logistic.results[,"AIC"]) 
    colnames(final_logistic)<-c("b1.ADD.OR","b1.P", "b2.DOMDEV OR","b2.P", "AIC") 
    rownames(final_logistic)<-paste("CHR",mapfile$V1,mapfile$V2,sep=" ")
    } else {
      final_logistic<-cbind(Logistic.results[,"b1.ADD"], b1.pvalues, Logistic.results[,"AIC"])
      colnames(final_logistic)<-c("b1.ADD.OR","b1.P", "AIC")
      rownames(final_logistic)<-paste("CHR",mapfile$V1,mapfile$V2,sep=" ")
  }
  #### Return the correct 2x3 or 2x2 table ###### Gen.table<-matrix()
  if (model == "genotypic") {
    Gen.table<-Table23Major(pedfile,mapfile) 
    } else if (model == "recessive") {
      tmp23<-Table23Major(pedfile,mapfile)$Genotypes 
      unaffecteds<-numeric(num.snps)
      for(i in 1:num.snps) {unaffecteds[i]<-tmp23[i,2]+tmp23[i,3]} 
      Gen.table<-cbind(tmp23[,1], unaffecteds) 
      colnames(Gen.table)<-c("AA", "Aa/aa")
    } else if (model == "dominant"){ 
      tmp23<-Table23Major(pedfile,mapfile)$Genotypes 
      affecteds<-numeric(num.snps)
      for(i in 1:num.snps) {
        affecteds[i]<-tmp23[i,1]+tmp23[i,2]} 
      Gen.table<-cbind(affecteds,tmp23[,3]) 
      colnames(Gen.table)<-c("AA/Aa","aa")
      } else { #additive/multiplicative Gen.table<-Table23Major(pedfile,mapfile)
      }
  return(list(REG=final_logistic, Table=Gen.table, LOG=logistic)) 
}
######################################################################
######################### END OF FUNCTIONS ###########################
######################################################################