#install.packages("multicore", repos="http://cran.utstat.utoronto.ca")

library(multicore)




reglev=function(x){
	y=y38[!is.na(x)]
	x1=x[!is.na(x)]

	chi2R=NA

	if(length(x1)>1){	
	if(length(table(x1))>2&min(table(x1))<=1){
		x1[x1==2]=1
	}
	if(length(table(x1))<3&min(table(x1))<=1){
		x1[x1==1]=0
		x1[x1==2]=0
	}
	}
	
	if(min(table(x1))>1){
		preg=summary(lm(y~x1))$coeff[2,4]
	}
	chi2R=qchisq(1-preg,1)
	return(cbind(chi2R)
}


reglevwrap<-function(xx){
	#y38=sample(new4$SAKNORM)
	y38=sample(new1$SaKnorm)
	assign("y38", y38, envir=.GlobalEnv)  # may or may not need this line
	return(apply(apply(new1[,seq(aa,bb)],2,reglev),1,sum))
}
#reglevwrap(1)

NN=10000

StatC=data.frame(do.call(rbind,mclapply(1:NN, reglevwrap,mc.cores=60)))

write.csv(StatC,'Observed_CAN_ALL-SLC9A3comp_Allvariants.csv',row.names=F)
StatC=read.csv('Observed_CAN_ALL-SLC9A3comp_Allvariants.csv')
head(StatC)

y38=(new1$SaKnorm)
chi2o=apply(apply(new1[,seq(aa,bb)],2,reglev),1,sum)

dim(StatC[StatC[,1]>chi2o[1],])[1]/dim(StatC)[1]
dim(StatC[StatC[,2]>chi2o[2],])[1]/dim(StatC)[1]
dim(StatC[StatC[,3]>chi2o[3],])[1]/dim(StatC)[1]




