library(mice)
library(scales)

load("../preprocess/proteome.young-old.n66.RData")
load("../resources/TF.CD.receptors.human.mouse.RData")

meanval<-rowMeans(data0,na.rm=T)
numval<-rowSums(!is.na(data))
numProt<-colSums(!is.na(data))
Q1<-quantile(meanval,seq(0,1,by=0.05))
ind5<- meanval>Q1["5%"]&meanval<Q1["95%"]

cum1<-cumsum(sort(numval,decreasing=T))

tempData <- mice(data[ind5&numval>39,],m=5,maxit=50,meth='pmm',seed=500)
completedData <- complete(tempData,1)

############# IMPUTATIONS ####################

completedDataY <- completedData[,indY] 
completedDataO <- completedData[,indO] 


res1<-prcomp(t(completedData))
res2<-prcomp((completedData))
res3<-prcomp(t(completedDataY ))
res4<-prcomp((completedDataY ))
res5<-prcomp(t(completedDataO ))
res6<-prcomp((completedDataO ))

pc1<-paste0(paste0("PC",seq(66)," ("),percent(res1$sdev^2/sum(res1$sdev^2)),")")
pc2<-paste0(paste0("PC",seq(66)," ("),percent(res2$sdev^2/sum(res2$sdev^2)),")")
pc3<-paste0(paste0("PC",seq(33)," ("),percent(res3$sdev^2/sum(res3$sdev^2)),")")
pc4<-paste0(paste0("PC",seq(33)," ("),percent(res4$sdev^2/sum(res4$sdev^2)),")")
pc5<-paste0(paste0("PC",seq(33)," ("),percent(res5$sdev^2/sum(res5$sdev^2)),")")
pc6<-paste0(paste0("PC",seq(33)," ("),percent(res6$sdev^2/sum(res6$sdev^2)),")")


COL_Compartments<-hue_pal()(4)[as.integer(factor(Compartments))]
PCH_ages<-c(17,16)[as.integer(factor(Ages))]
COL_levels<-seq(3)[as.integer(factor(Levels))]
COL_IO<-hue_pal()(2)[as.integer(factor(innerOuter))]


LPCA<-list(res1,res2,res3,res4,res5,res6)
Lpercent<-list(pc1,pc2,pc3,pc4,pc5,pc6)

########################################

source("step2.help.funcs.R")
pdf("PCA.66.with.SVM.DEC17.pdf")
	Q2<-quantile(meanval)
	hist(meanval,breaks=100,main="histogram of average values per gene\nfor all 3100 genes")
	abline(v=quantile(meanval),col=4,lwd=2,lty=2)
	text(Q2,sample(seq(20,120),length(Q2)),names(Q2),srt=90)

	hist(meanval,breaks=100,main="histogram of average values per gene\nfor all 3100 genes")
	Q1<-quantile(meanval,seq(0,1,by=0.05))
	abline(v=Q1,col=2,lwd=2,lty=2)
	abline(v=quantile(meanval),col=4,lwd=2,lty=2)
	text(Q1,sample(seq(20,120),length(Q1)),names(Q1),srt=90)


	par(mfrow=c(2,1))
	y<-rev(sort(numval))
	z<-c(rep(0,200),y[-seq(200)]-y[seq(length(y)-200)])

	bpp<-barplot(y,ylab="Number of valid values",xlab="genes (in decreasing order of # of valid values)")
	PEAK<-max(which(z==min(z)))
	abline(v=bpp[PEAK],lwd=2,col=2,lty=2)
	abline(h=y[PEAK],lwd=2,col=2,lty=2)
	text(bpp[PEAK],66,PEAK,srt=0,col=2,xpd=T)
	legend("topright",inset=c(0,-0.05,0,0),
		legend="the optimal cutoff is '# valid values>39 samples'.
			it corresponds to the steepest slope. 
			beyond this point, miss values are increasing 
			faster than valid values.",cex=2/3,xpd=T,bty='n')

	barplot(z,ylab="slope")
	abline(v=bpp[PEAK],lwd=2,col=2,lty=2)
	text(bpp[PEAK],0,PEAK,srt=0,col=2,xpd=T)


	x<-cum1/max(cum1)*100

	r1<-100-cum1/(seq(length(cum1))*66)*100

	bpp<-barplot(x,ylab="% of all valid values",main="% of all valid values in whole dataset")
	abline(v=bpp[PEAK],lwd=2,col=2,lty=2)
	abline(h=c(x[PEAK],100),lwd=2,col=2,lty=2)
	text(bpp[PEAK]*1.5,x[PEAK]*1.2,percent(x[PEAK]/100),srt=0,col=2,xpd=T)

	bpp<-barplot(r1,ylab="% missing value if cutoff is made here",
		main="% of missing values for each cutoff")
	abline(v=bpp[PEAK],lwd=2,col=2,lty=2)
	abline(h=c(r1[PEAK],100),lwd=2,col=2,lty=2)
	text(bpp[PEAK]*1.5,r1[PEAK]*1.2,percent(r1[PEAK]/100),srt=0,col=2,xpd=T)

	z<-c(rep(0,200),x[-seq(200)]-x[seq(length(x)-200)])
	PEAK<-max(which(z==max(z)))
	#barplot(z,ylab="slope")
	#abline(v=bpp[PEAK],lwd=2,col=2,lty=2)
	#text(bpp[PEAK],0,PEAK,srt=0,col=2,xpd=T)

	par(mfrow=c(1,1))



	plot(numval,meanval,xlab="num of valid values, per gene",
		ylab="average expr levels among valid values, per gene")
	abline(h=mean(meanval),lty=2,col=4)
	points(numval[ind5],meanval[ind5],pch=16,col="red")
	legend("topleft",pch=c(16,1),cex=1,col=c(2,1),
		legend=c(paste0("avrg expr in the 5%~95% range (",sum(ind5)," genes)"),
		paste0("avrg expr outside 5%~95% range (",sum(!ind5)," genes)")))
###
	source("step2.help.funcs.R")

	ind1<-ind5&numval>39
	PCAfunc(ind1,"(N/A imputed)",1,"39 samples& in 5%~95% range",indS=seq(66),T,T,T)

	ind1<-ind5&numval>39
	PCAfunc(ind1,"(N/A imputed)",2,"39 samples& in 5%~95% range",indY,T,T,T)

	ind1<-ind5&numval>39
	PCAfunc(ind1,"(N/A imputed)",3,"39 samples& in 5%~95% range",indO,T,T,T)

dev.off()


