setwd("/Users/naim panjwani/Documents/Strug/UKRE_assoc/")

logistic <- read.table("02_sorted_logistic.logistic",header=T)

U<-sort(runif(length(logistic$P)))

qqplot(-log(logistic$P,base=10),-log(U,base=10),xlab="-log10 P-values",ylab="-log10 Expected P-values", main="Q-Q Plot Logistic Regression\nAdditive Model")
lines(c(0,5),c(0,5))

