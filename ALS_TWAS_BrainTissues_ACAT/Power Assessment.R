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
library(MASS)
corv = as.matrix(cor1)
diag(corv) = 1

n=1e3
tx = seq(1e-5,0.05,length=20)
set.seed(1234567)
mu=matrix(rnorm(13*n,0,2.5),n,13)
zx=mu
for (i in 1:n){zx[i,] = mvrnorm(1,mu=mu[i,],Sigma=corv)}

pvalue = function(x) {-pnorm(abs(x))*2}
ax = zx
for (j in 1:13)
{
ax[,j] = pnorm(-abs(zx[,j]))*2
}

fdr = function(x) {p.adjust(x, method ="bonferroni")}
res1=ax
p=dim(res1)[2]
px=matrix(NA,length(tx),p)
for (j in 1:length(tx))
{
for (i in 1:dim(res1)[2])
{
px[j,i]= mean((na.omit(res1[,i])) < (tx[j]/dim(res1)[1]))
}
}
matplot(tx, px, type="l",lwd=2,col=c(1:p)+1,ylim=c(0,0.85),pch=15,ylab="Estimated power",xlab="",main="")
for (j in 1:p)
{
points(tx,px[,j],pch=15,col=j+1,cex=0.1)
}
for (j in 1:1)
{
lines(tx,px[,j],pch=15,lwd=3,col=j+1,cex=0.1)
}
for (i in 1:5){mtext(expression(paste("(E)",sep="")),outer=F,at=-0.01,cex=1.2)}
if (S==5) for (i in 1:5){legend("topleft",c("ACAT13","ACAT8","ACAT4","ACAT1","oracle"),
	bty="n",pch=c(15:19),col=c(1:(p+1))+1,cex=1)}



CX0 = solve(corv)
n=1e3
tx = seq(1e-5,0.05,length=20)
set.seed(1234567)
mu=matrix(rnorm(13*n,0,2.5),n,13)
zx=mu
for (i in 1:n){zx[i,] = mvrnorm(1,mu=mu[i,],Sigma=corv)}
pvalue = function(x) {-pnorm(abs(x))*2}
ax = zx
for (j in 1:13)
{
ax[,j] = pnorm(-abs(zx[,j]))*2
}

p=dim(zx)[1]
ux=rep(NA,p)
cx = 3
for (j in 1:p)
{
idm = sort(sample(seq(1:dim(zx)[2]),cx,replace=F))
z = t(zx[j,idm]) %*% CX0[idm,idm] %*% (zx[j,idm])
ux[j] = pchisq(z,13,lower.tail=FALSE)
}


res=ax
p=dim(res)[1]
acat1=rep(NA,p)
acat2=rep(NA,p)
acat3=rep(NA,p)
acat4=rep(NA,p)
for (j in 1:p)
{
pxx=c(na.omit(unlist(res[j,])))
acat1[j] = ACAT(pxx)
idm = sort(sample(seq(1:dim(res)[2]),4,replace=F))
pxx = res[j,]
pxx[idm] = NA
pxx=c(na.omit(unlist(pxx)))
acat2[j] = ACAT(pxx)
idm = sort(sample(seq(1:dim(res)[2]),8,replace=F))
pxx = res[j,]
pxx[idm] = NA
pxx=c(na.omit(unlist(pxx)))
acat3[j] = ACAT(pxx)
idm = sort(sample(seq(1:dim(res)[2]),12,replace=F))
pxx = res[j,]
pxx[idm] = NA
pxx=c(na.omit(unlist(pxx)))
acat4[j] = ACAT(pxx)
}

fdr = function(x) {p.adjust(x, method ="bonferroni")}
res1=cbind(acat1,acat2,acat3,acat4,ux)
p=dim(res1)[2]
px=matrix(NA,length(tx),p)
for (j in 1:length(tx))
{
for (i in 1:dim(res1)[2])
{
px[j,i]= mean((na.omit(res1[,i])) < (tx[j]/dim(res1)[1]))
}
}
matplot(tx, px, type="l",lwd=2,col=c(1:p)+1,ylim=c(0,0.85),pch=15,ylab="Estimated power",xlab="",main="")
for (j in 1:p)
{
points(tx,px[,j],pch=15,col=j+1,cex=0.1)
}
for (j in 1:1)
{
lines(tx,px[,j],pch=15,lwd=3,col=j+1,cex=0.1)
}
for (i in 1:5){mtext(expression(paste("(F)",sep="")),outer=F,at=-0.01,cex=1.2)}
