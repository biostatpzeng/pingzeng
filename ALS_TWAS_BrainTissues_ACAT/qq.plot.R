# from the internet
qq = function(pvector,conf.points=10000,colx,conf.alpha=.05, maxx,...) {
	pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
	o = -log10(sort(pvector,decreasing=F))
	e = -log10(ppoints(length(pvector)))
	if (is.null(maxx)) maxx = min(max(o,e),8)*1.05
	plot(c(0,100),c(0,100),xlim=c(0,maxx),ylim=c(0,maxx),pch=1,cex=0,
		xlab=expression(Expected~~-log[10](italic(p))),
		ylab=expression(Observed~~-log[10](italic(p))),
		)
	n <- length(pvector)+1
	conf.points = min(conf.points, n-1);
	mpts<-matrix(nrow=conf.points*2, ncol=3)
	for(i in seq(from=1, to=conf.points)) {
	mpts[i,1]<- -log10((i-.5)/n)
	mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
	mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
	mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
	}
	polygon(x=mpts[,1],y=mpts[,2], col="gray",lty=1,lwd =0.1,border= "gray")
	points(e,o,pch=20,cex=0.5,col=colx)
	abline(0,1,col="red",lty=2)
}

