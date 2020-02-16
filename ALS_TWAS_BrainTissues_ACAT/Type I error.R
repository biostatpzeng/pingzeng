library(data.table)
tx=c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24",
"Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9",
"Brain_Hippocampus","Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
"Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")
genex = NULL
for (j in 1:length(tx))
{
ax=data.frame(fread(paste(tx[j],".txt",sep=""),sep='\t',head=TRUE))
genex=c(genex,ax$ID)
}

library(data.table)
tx=c("Brain_Amygdala","Brain_Anterior_cingulate_cortex_BA24",
"Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9",
"Brain_Hippocampus","Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
"Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra")
genex1=unique(genex)
px = matrix(NA,length(genex1),13)
for (j in 1:length(tx))
{
ax=data.frame(fread(paste(tx[j],".txt",sep=""),sep='\t',head=TRUE))
m = length(ax$ID)
index = rep(NA,m)
for (i in 1:m)
{
index[i] = which(genex1==ax$ID[i])
}
px[index,j]=ax$TWAS.Z
print(j)
}

cor1=matrix(NA,13,13)
for (i in 1:13)
{
for (j in 1:13)
{
cor1[i,j] = cor(na.omit((px[,c(i,j)])))[1,2]
}
}

rownames(cor1)=c("Amygdala","Anterior cingulate cortex BA24",
"Caudate basal ganglia","Cerebellar Hemisphere","Cerebellum","Cortex",
"Frontal Cortex BA9","Hippocampus","Hypothalamus","Nucleus accumbens basal ganglia",
"Putamen basal ganglia","Spinal cord cervical","Substantia nigra")


source("qq.plot.R")
source("ACAT_function.R")
library(harmonicmeanp)
library(MASS)
corv = as.matrix(cor1)
corv = as.matrix(cor1)*0
diag(corv) = 1

n=1e6
zx = mvrnorm(n,mu=rep(0,13),Sigma=corv,empirical=TRUE)
pvalue = function(x) {-pnorm(abs(x))*2}
px = zx
for (j in 1:13)
{
px[,j] = pnorm(-abs(zx[,j]))*2
}

res=px
p=dim(res)[1]
acat=rep(NA,p)
fihser=rep(NA,p)
for (j in 1:p)
{
pxx=c(na.omit(unlist(res[j,])))
acat[j] = ACAT(pxx)
fihser[j] = pchisq((-sum(log(pxx))*2), df=length(pxx)*2, lower.tail=F)
}

par(mfrow=c(2,2), cex.axis=1.5, cex.lab=1.5, cex.main=1.5, lwd=1,
    mgp=c(2.5, 0.5, 0), tcl=-0.2, font.axis=1.5, font.lab=1.5, 
    mar=c(5.1,6,2.5,2.1))

qq(acat,main="",cex=0.1,pch=15,colx="red")
pvector=fihser
o = -log10(sort(pvector,decreasing=F))
e = -log10(ppoints(length(pvector)))
points(e,o,pch=16,cex=0.1,col="blue")
for (i in 1:5){mtext(expression(paste("(A)",sep="")),outer=F,at=-0.5,cex=1.2)}


res=px
p=dim(res)[1]
acat=rep(NA,p)
fihser=rep(NA,p)
for (j in 1:p)
{
idm = sort(sample(seq(1:dim(res)[2]),round(runif(1,5,dim(res)[2]-2)),replace=F))
res[j,idm] = NA
pxx=c(na.omit(unlist(res[j,])))
acat[j] = ACAT(pxx)
fihser[j] = pchisq((-sum(log(pxx))*2), df=length(pxx)*2, lower.tail=F)
}

qq(acat,main="",cex=0.1,pch=15,colx="red")
pvector=fihser
o = -log10(sort(pvector,decreasing=F))
e = -log10(ppoints(length(pvector)))
points(e,o,pch=16,cex=0.1,col="blue")
for (i in 1:5){mtext(expression(paste("(B)",sep="")),outer=F,at=-0.5,cex=1.2)}



source("ACAT_function.R")
library(MASS)
corv = as.matrix(cor1)
diag(corv) = 1
n=1e6
zx = mvrnorm(n,mu=rep(0,13),Sigma=corv,empirical=TRUE)
pvalue = function(x) {-pnorm(abs(x))*2}
px = zx
for (j in 1:13)
{
px[,j] = pnorm(-abs(zx[,j]))*2
}

res=px
p=dim(res)[1]
acat=rep(NA,p)
fihser=rep(NA,p)
for (j in 1:p)
{
pxx=c(na.omit(unlist(res[j,])))
acat[j] = ACAT(pxx)
fihser[j] = pchisq((-sum(log(pxx))*2), df=length(pxx)*2, lower.tail=F)
}

qq(acat,main="",cex=0.1,pch=15,colx="red")
pvector=fihser
o = -log10(sort(pvector,decreasing=F))
e = -log10(ppoints(length(pvector)))
points(e,o,pch=16,cex=0.1,col="blue")
for (i in 1:5){mtext(expression(paste("(C)",sep="")),outer=F,at=-0.5,cex=1.2)}

res=px
p=dim(res)[1]
acat=rep(NA,p)
fihser=rep(NA,p)
for (j in 1:p)
{
idm = sort(sample(seq(1:dim(res)[2]),round(runif(1,5,dim(res)[2]-2)),replace=F))
res[j,idm] = NA
pxx=c(na.omit(unlist(res[j,])))
acat[j] = ACAT(pxx)
fihser[j] = pchisq((-sum(log(pxx))*2), df=length(pxx)*2, lower.tail=F)
}

qq(acat,main="",cex=0.1,pch=15,colx="red")
pvector=fihser
o = -log10(sort(pvector,decreasing=F))
e = -log10(ppoints(length(pvector)))
points(e,o,pch=16,cex=0.1,col="blue")
for (i in 1:5){mtext(expression(paste("(D)",sep="")),outer=F,at=-0.5,cex=1.2)}
for (i in 1:5){legend("bottomright",c("Fisher's","ACAT"),bty="n",pch=c(15),col=c("blue","red"),cex=1.5)}



