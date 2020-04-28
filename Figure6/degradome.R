library(scales)

NP<-read.csv("NP-degr.txt",sep="\t")
AF<-read.csv("AF-degr.txt",sep="\t")

###############
indCDSNP<-grep("^CDS.*NP",colnames(NP))
indDSNP<-grep("^DS.*NP",colnames(NP))
indmedNP<-grep("^medi",colnames(NP))

CDSNP<-NP[,indCDSNP[-1]]
DSNP<-NP[,indDSNP]

CDSNPsd<-apply(CDSNP,1,sd,na.rm=T)
CDSNPsd[is.na(CDSNPsd)]<-0

DSNPsd<-apply(DSNP,1,sd,na.rm=T)
DSNPsd[is.na(DSNPsd)]<-0
SEQNP<-seq(nrow(CDSNP))

###############
indCDSAF<-grep("^CDS.*AF",colnames(AF))
indDSAF<-grep("^DS.*AF",colnames(AF))
indodd<-seq(1,nrow(AF),by=2)
ODDAF<-AF[indodd,]

indmedAF<-grep("^medi",colnames(AF))
ODDAF[,indmedAF]

CDSAF<-ODDAF[,indCDSAF[-1]]
DSAF<-ODDAF[,indDSAF]

CDSAFsd<-apply(CDSAF,1,sd,na.rm=T)
CDSAFsd[is.na(CDSAFsd)]<-0

DSAFsd<-apply(DSAF,1,sd,na.rm=T)
DSAFsd[is.na(DSAFsd)]<-0
SEQAF<-seq(nrow(CDSAF))

###############
pdf("degradome.apr26.pdf",height=4,width=12)
	layout(t(matrix(seq(2))),widths=c(9,3))
	bpp<-matplot(cbind(CDSNP,DSNP),pch=16,type='n',
		cex=2,col=rep(rev(hue_pal()(2)),c(2,3)),
		axes=F)
	axis(2)
	lines(SEQNP,NP[,indmedNP],lwd=2)
	polygon(c(SEQNP,rev(SEQNP)),
		c(NP[,"CDS.average.normalised"]+CDSNPsd,
		rev(NP[,"CDS.average.normalised"]-CDSNPsd)),
		border=NA,col="#00BFC460")
	lines(SEQNP,NP[,"CDS.average.normalised"],col="#00BFC4",lwd=3)

	polygon(c(SEQNP,rev(SEQNP)),
		c(NP[,"DS.average.normalised"]+DSNPsd,
		rev(NP[,"DS.average.normalised"]-DSNPsd)),
		border=NA,col="#F8766D60")
	lines(SEQNP,NP[,"DS.average.normalised"],col="#F8766D",lwd=3)
	bpp<-matplot(cbind(CDSNP,DSNP),pch=16,
		cex=1,col=rep(rev(hue_pal()(2)),c(2,3)),add=T)
	abline(h=0)
	axis(1,at=SEQNP,tick=T,
		labels=gsub(".*GN=| .*","",NP[,2]),cex.axis=1/2,
		las=2)
	#x<-as.vector(as.matrix(CDSNP))
	#y<-as.vector(as.matrix(DSNP))
	x<-NP[,"CDS.average.normalised"]
	y<-NP[,"DS.average.normalised"]

 t.test(y-x)
	indy<-which(y>x)
	fitx <- density(x[indy])
	fity <- density(y[indy])
	plot(fitx, xlim=c(-2,6))
	lines(fity$x,fity $y,col=2,lwd=2)


 pie(table(sign(y-x)))

	fitx <- density(x)
	fity <- density(y)
	plot(fitx,xlim=c(-6,4))
	lines(fity$x,fity $y,col=2,lwd=2)
library(ROCit)
## Warning: package 'ROCit' was built under R version 3.5.2
DAT<-data.frame(pred=c(x,y),Y=rep(c("YND","AGD"),c(length(x),length(y))))
ROCit_obj <- rocit(score=DAT$pred,class=DAT$Y)
res1<-plot(ROCit_obj)

###########
	bpp<-matplot(cbind(CDSAF,DSAF),pch=16,type='n',
		cex=2,col=rep(rev(hue_pal()(2)),c(2,3)),
		axes=F)
	axis(2)
	lines(SEQAF,ODDAF[,indmedAF],lwd=2)
	polygon(c(SEQAF,rev(SEQAF)),
		c(ODDAF[,"CDS.AVERAGE.NORMAL"]+CDSAFsd,
		rev(ODDAF[,"CDS.AVERAGE.NORMAL"]-CDSAFsd)),
		border=NA,col="#00BFC460")
	lines(SEQAF,ODDAF[,"CDS.AVERAGE.NORMAL"],col="#00BFC4",lwd=3)

	polygon(c(SEQAF,rev(SEQAF)),
		c(ODDAF[,"DS.AVERAGE.NORMAL"]+DSAFsd,
		rev(ODDAF[,"DS.AVERAGE.NORMAL"]-DSAFsd)),
		border=NA,col="#F8766D60")
	lines(SEQAF,ODDAF[,"DS.AVERAGE.NORMAL"],col="#F8766D",lwd=3)
	bpp<-matplot(cbind(CDSAF,DSAF),pch=16,
		cex=1,col=rep(rev(hue_pal()(2)),c(2,3)),add=T)
	abline(h=0)
	axis(1,at=SEQAF,tick=T,
		labels=gsub(".*GN=| .*","",ODDAF[,2]),cex.axis=1/2,
		las=2)
	#x<-as.vector(as.matrix(CDSAF))
	#y<-as.vector(as.matrix(DSAF))
	#y<-y[!is.na(y)]
	x<-ODDAF[,"CDS.AVERAGE.NORMAL"]
	y<-ODDAF[,"DS.AVERAGE.NORMAL"]
 pie(table(sign(y-x)))

	indy<-which(y>x)
	fitx <- density(x[indy])
	fity <- density(y[indy])
	plot(fitx, xlim=c(-2,6))
	lines(fity$x,fity $y,col=2,lwd=2)
dev.off()


