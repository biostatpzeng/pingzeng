
	library(data.table)
	library(MendelianRandomization)

#---------------------------------------------------Main analysis-------------------------------------------------#
	yx = read.table("44indexdrink_ALS.txt",header = T,sep = "\t")
	MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
	IVWObject = mr_ivw(MRinput,model="fixed")

#	> IVWObject
#	Inverse-variance weighted method
#	(variants uncorrelated, fixed-effect model)
#
#	Number of Variants : 44 
#
#	------------------------------------------------------------------
#	 Method Estimate Std Error 95% CI       p-value
#	    IVW    0.907     0.298 0.324, 1.491   0.002
#	------------------------------------------------------------------
#	Residual standard error =  1.137 
#	Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.
#	Heterogeneity test statistic = 55.5512 on 43 degrees of freedom, (p-value = 0.0950)

	EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)

#	> EggerObject
#
#	MR-Egger method
#	(variants uncorrelated, random-effect model)
#
#	Number of Variants =  44 
#
#	------------------------------------------------------------------
#	      Method Estimate Std Error  95% CI       p-value
#	    MR-Egger    0.784     1.695 -2.638, 4.205   0.646
#	 (intercept)    0.001     0.013 -0.026, 0.028   0.941
#	------------------------------------------------------------------

	weighted_MedianObject <- mr_median(MRinput, weighting = "weighted",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)

#	> weighted_MedianObject

#	 Weighted median method 
#
#	Number of Variants : 44 
#	------------------------------------------------------------------
#			 Method Estimate Std Error  95% CI       p-value
#	 Weighted median method    0.678     0.445 -0.194, 1.551   0.127
#	------------------------------------------------------------------


#------------------------------------------------------loocv------------------------------------------------------#
	yx = read.table("44indexdrink_ALS.txt",header = T,sep = "\t")	
	MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
	IVWObject = mr_ivw(MRinput,model="fixed")
	EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)

	res = matrix(NA,dim(yx)[1],4)
	for( z in 1:dim(yx)[1]){
	  temp = yx[-z,]
	  MRinput1 = mr_input(bx = temp$BETA.x ,bxse = temp$SE.x,by = temp$BETA.y,byse = temp$SE.y)
	  IVWObject1 = mr_ivw(MRinput1,model="fixed")
	  res[z,] = c(IVWObject1$Estimate,1/IVWObject1$CILower,1/IVWObject1$CIUpper,IVWObject1$Pvalue)
	}
	res = data.frame(yx$SNP,res)
	colnames(res) = c("SNP","BETA","LOW","UP","PVAL")


#-------------------------------------------------remove "rs112635299"--------------------------------------------#
	yx = read.table("44indexdrink_ALS.txt",header = T,sep = "\t")
	index = which(yx$SNP%in%c("rs112635299"));if (length(index)>0){yx=yx[-index,]}
	MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
	IVWObject = mr_ivw(MRinput,model="fixed")
	EggerObject <- mr_egger(MRinput,robust = FALSE,penalized = FALSE,correl = FALSE,distribution = "t-dist",alpha = 0.05)
	MaxlikObject <- mr_maxlik(MRinput,correl = FALSE,distribution = "normal",alpha = 0.05)
	weighted_MedianObject <- mr_median(MRinput, weighting = "weighted",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)
	simple_MedianObject <- mr_median(MRinput, weighting = "simple",distribution = "normal", alpha = 0.05, iterations = 10000,seed = 314159265)


#-------------------------------------------------remove EBI SNPs-------------------------------------------------#
	yx = read.table("44indexdrink_ALS.txt",header = T,sep = "\t")
	index = which(yx$SNP%in%c("rs13107325","rs9320010"));if (length(index)>0){yx=yx[-index,]}
	ebi_snp = as.character(read.table("EBI_SNP.txt",header = F)[,1])
	index = which(yx$SNP%in%ebi_snp)
	yx = yx[-index,]
	MRinput = mr_input(bx = yx$BETA.x ,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y)
	IVWObject = mr_ivw(MRinput,model="fixed")

#	> IVWObject
#
#	Inverse-variance weighted method
#	(variants uncorrelated, fixed-effect model)
#
#	Number of Variants : 29 
#
#	------------------------------------------------------------------
#	 Method Estimate Std Error 95% CI       p-value
#	    IVW    0.823     0.381 0.076, 1.570   0.031
#	------------------------------------------------------------------
#	Residual standard error =  1.203 
#	Residual standard error is set to 1 in calculation of confidence interval by fixed-effect assumption.
#	Heterogeneity test statistic = 40.5398 on 28 degrees of freedom, (p-value = 0.0591)


#-----------------------------------------------multivariable MR--------------------------------------------------#

	yx1 = read.table("index_snp_multivariable.txt",header = T,sep = "\t")
	res = summary(lm(yx1$BETA.y ~ yx1$BETA.x + yx1$BETA - 1,weights = 1/(yx1$SE.y*yx1$SE.y)))
#	Coefficients:
#		   Estimate Std. Error t value Pr(>|t|)  
#	yx1$BETA.x   0.8030     0.3794   2.116   0.0403 *
#	yx1$BETA    -0.3275     0.5228  -0.626   0.5344  


#-----------------------------------------------UKB types of alcohol----------------------------------------------#

	x1=c("Freq","Red","ChWhite","BeCider","spirits","Fortified","WithM","10Years","former","never","previous")
	m = length(x1)
	res=matrix(NA,m,4)
	for (j in 1:m){
		try({
			data = read.table(paste("Additional_alcohol_index_snp.txt",sep=""),head=T,fill=T)
			index_type = which(data$TYPE == x1[j])
			yx = data[index_type,]
			index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
			if (length(index1)<length(yx$BETA.x)){yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
			index2 = which(yx$P.y <(0.05/dim(yx)[1]))
			if (length(index2)>0) {yx=yx[-index2,]}
			yx=yx[which(yx$P.x<1E-6),]
			fit = mr_ivw(mr_input(bx = yx$BETA.x,bxse = yx$SE.x,by = yx$BETA.y,byse = yx$SE.y),model = "fixed")
			res[j,]=c(dim(yx)[1],fit$Estimate,fit$StdError,fit$Pvalue)
			print(j)
		})
		}


