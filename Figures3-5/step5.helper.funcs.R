source("../resources/CAT_DEGs.R")

load("../resources/MATRISOME.RData")
load("../resources/NABA.break.down.RData")
load("../resources/TF.CD.receptors.human.mouse.RData")

LNABA[["Collagens"]]<-as.character(MATRISOME[grep("^COL[0-9]",MATRISOME[,1]),1])
LNABA2<-LNABA
LNABA2[["TF"]]<-TF.human
LNABA2[["CD"]]<-surface.human

library(grid)
library(lattice)
library(gridExtra)
library(stringr)

load("../resources/anita.RData")
indAnita<-which(geneSymb2%in%gsub(" .*","",annots[,2]))
indAnita2<-match(geneSymb2,gsub(" .*","",annots[,2]))
#indAnita2<-match(intersect(geneSymb2,MATRISOME[,1]),gsub(" .*","",annots[,2]))
ALOG<-dat.log[indAnita2,]
yAF<-ALOG[,c(1,3)]
yNP<-ALOG[,c(2,4)]
oAF<-ALOG[,c(5,7)]
oNP<-ALOG[,c(6,8)]
yAFNPdiff<-rowMeans(yAF)-rowMeans(yNP)
oAFNPdiff<-rowMeans(oAF)-rowMeans(oNP)
oyAFdiff<-rowMeans(oAF)-rowMeans(yAF)
oyNPdiff<-rowMeans(oNP)-rowMeans(yNP)
Ldiff<-list("young AF - young NP (microarray)"=yAFNPdiff,
	"old AF - old NP (microarray)"=oAFNPdiff,
	"old AF - young AF (microarray)"=oyAFdiff,
	"old NP - young NP (microarray)"=oyNPdiff)
withANITA<-function(xvar,res1){
	par(mfrow=c(1,2))
	cond1<-res1$cond1
	cond2<-res1$cond2
	indS1<-res1$indS1
	indS2<-res1$indS2
	indT<-res1$indT
	pMAT<-res1$pMAT
	indSig<-res1$indSig
	indSig2<-res1$indSig2
	DIFF<-pMAT[2,indT][indSig]

	print(str(Ldiff[[xvar]]))
	m1<-rowMeans(data[,indS1],na.rm=T)
	m2<-rowMeans(data[,indS2],na.rm=T)
	plot(Ldiff[[xvar]],m2-m1,main=paste0("comparion ",compi,": Compared with Anita"),
		pch=".",xlab=xvar,ylab=paste0(cond2," - ",cond1),
		col="grey")
	abline(h=0,col=3,lty=2)
	abline(v=0,col=3,lty=2)
	if(length(indSig)>0){
		points(Ldiff[[xvar]][indT][indSig][DIFF>0],(m2-m1)[indT][indSig][DIFF>0],
			pch=16,xlab=xvar,ylab=paste0(cond2," - ",cond1),
			cex=2,col="#F8766D")
		points(Ldiff[[xvar]][indT][indSig][DIFF<0],(m2-m1)[indT][indSig][DIFF<0],
			pch=16,xlab=xvar,ylab=paste0(cond2," - ",cond1),
			cex=2,col="#619CFF")
		indshow1<- (Ldiff[[xvar]][indT][indSig]>0 & (m2-m1)[indT][indSig]>0)
		indshow2<- (Ldiff[[xvar]][indT][indSig]<0 & (m2-m1)[indT][indSig]<0)
		indshow12<- indshow1|indshow2
		tab1<-table(Ldiff[[xvar]][indT][indSig]>0,(m2-m1)[indT][indSig]>0)
		tab2<-table(Ldiff[[xvar]][indT][indSig]<0,(m2-m1)[indT][indSig]<0)
		print(tab1)
		print(tab2)
		text(Ldiff[[xvar]][indT][indSig][indshow12],
			(m2-m1)[indT][indSig][indshow12],
			geneSymb2[indT][indSig][indshow12],cex=1/2)
	}
	par(mfrow=c(1,2))
}

DRAWGENE<-function(downDEG,upDEG){
		L1<-sapply(LNABA2,function(x)intersect(x,downDEG))
		L1[["others"]]<-setdiff(downDEG,unlist(L1))
		str1<-sapply(L1,paste0,collapse=", ")

		L2<-sapply(LNABA2,function(x)intersect(x,upDEG))
		L2[["others"]]<-setdiff(upDEG,unlist(L2))
		str2<-sapply(L2,paste0,collapse=", ")

		mat1<-cbind(str_wrap(str1,40),str_wrap(str2,40))
		colnames(mat1)<-c(paste0("Down in RED (Up in Blue);n=",length(downDEG)),
			paste0("Up in RED (Down in Blue);n=",length(upDEG)))
		mytheme <- gridExtra::ttheme_default(base_size =6,padding = unit(c(1, 1), "mm"))
		rownames(mat1)<-paste0(names(L1)," (",sapply(L1,length),";",sapply(L2,length),")")

		pushViewport(viewport(x=0.75,y=0.5))
		grid.table(mat1,theme=mytheme )
		popViewport()
}
ONOFF<-function(vecIN,cond1="young",cond2="old"){
	indS1<-which(vecIN%in%cond1)
	indS2<-which(vecIN%in%cond2)
	m1<-rowMeans(data[,indS1],na.rm=T)
	m2<-rowMeans(data[,indS2],na.rm=T)

	flag1<-rowSums(!is.na(data[,indS1]))
	flag2<-rowSums(!is.na(data[,indS2]))

	indUp  <-which(flag1>1 & m1>0 & flag2==0)
	indDown<-which(flag2>1 & m2>0 & flag1==0)

	indUp2  <-which(flag1>(1/2*length(indS1)) & m1>0 & flag2==0)
	indDown2<-which(flag2>(1/2*length(indS2)) & m2>0 & flag1==0)


	indSig<-c(indUp,indDown)
	M1<-c(m1[indUp],-m2[indDown])
	numY<-c(flag1[indUp],flag2[indDown])

	XLAB1<-paste0("# samples on in [",cond1,"]")
	XLAB2<-paste0("# samples on in [",cond2,"]")

	YLAB1<-paste0("avrg expression in [",cond1,"]")
	YLAB2<-paste0("avrg expression in [",cond2,"]")
	MAIN1<-paste0("Comparison ",compi,"-D: non-statistical DEGs\nturned on one-sided in [",cond1,"]")
	MAIN2<-paste0("Comparison ",compi,"-D: non-statistical DEGs\nturned on one-sided in [",cond2,"]")

	par(mfrow=c(1,3))
	layout(t(matrix(seq(3))),widths=c(4,4,4))

	if(length(indUp)>0){
		plot(flag1[indUp],m1[indUp],xlab=XLAB1,ylab=YLAB1,main=MAIN1,pch=".",col="#619CFF",xlim=c(0,length(indS1)))
		abline(v=1/2*length(indS1),lth=2)
		if(length(indUp2)>0)
			text(flag1[indUp2],m1[indUp2],geneSymb2[indUp2],xpd=T,cex=1/2)
	}

	if(length(indDown)>0){
		plot(flag2[indDown],m2[indDown],xlab=XLAB2,ylab=YLAB2,main=MAIN2,pch=".",col="#F8766D",xlim=c(0,length(indS2)))
		abline(v=1/2*length(indS2),lth=2)
		if(length(indDown2)>0)
			text(flag2[indDown2],m2[indDown2],geneSymb2[indDown2],xpd=T,cex=1/2)
	}
	BSIZE<- 5 -length(c(indDown,indUp))/3e2
	if(BSIZE<3)BSIZE<-3

	CAT_DEG2(geneSymb2[indDown],geneSymb2[indUp],BSIZE=4,Xpos=0.82,WIDTH=60)
	if(exists("OUTDIR")){
		COND1<-paste0(cond1,".n",length(indS1))
		COND2<-paste0(cond2,".n",length(indS2))
		if(length(c(indUp2,indDown2))>0){
			datout<-data.frame(gene=c(geneSymb2[indUp2],geneSymb2[indDown2]),
				status=rep(c(COND1,COND2),c(length(indUp2),length(indDown2))),
				COND1=COND1,
				COND2=COND2,
				len1=length(indS1),
				len2=length(indS2),
				numON1=c(flag1[indUp2],flag1[indDown2]),
				numON2=c(flag2[indUp2],flag2[indDown2]))
			prefix<-paste0("No.",formatC(compi,width=3,flag = "0"),".",COND1,"-vs-",COND2,".nDEG",nrow(datout))
			fnout=paste0(OUTDIR,"/One-sied.",prefix,".txt")
			cat("writing one-sided results to: ",fnout,"\n")
			print(head(datout))
			write.table(datout,file=fnout,sep="\t",quote=F,row.names=F)
		}
	}

}
FOLD_CHNAGE<-function(vecIN,cond1="young",cond2="old",indT,indSig){
	indS1<-which(vecIN%in%cond1)
	indS2<-which(vecIN%in%cond2)
	m1<-rowMeans(data[,indS1],na.rm=T)
	m2<-rowMeans(data[,indS2],na.rm=T)

	flag1<-rowSums(!is.na(data[,indS1]))
	flag2<-rowSums(!is.na(data[,indS2]))

	indFC1<-(m1-m2)>2 & (m1+m2)>20
	indFC2<-(m1-m2)< -2 & (m1+m2)>20

	par(mfrow=c(1,2))
	plot(m1,m2,pch=".",col="grey",
		main=paste0("Comparison:",compi,"-C: non-statistical DEGs based on FC"),
		xlab=cond1,ylab=cond2)
	if(length(indFC1)>0)points(m1[indFC1],m2[indFC1],pch=16,col="#619CFF")
	if(length(indFC2)>0)points(m1[indFC2],m2[indFC2],pch=16,col="#F8766D")

	indFC21<-setdiff(which(indFC1),which(indT)[indSig])
	indFC22<-setdiff(which(indFC2),which(indT)[indSig])
	abline(0,1)

	if(length(indFC21)>0)text(m1[indFC21],m2[indFC21],geneSymb2[indFC21],cex=1/2)
	if(length(indFC22)>0)text(m1[indFC22],m2[indFC22],geneSymb2[indFC22],cex=1/2)
	CAT_DEG2(geneSymb2[indFC22],geneSymb2[indFC21], BSIZE=5, WIDTH=60,Xpos=0.75)
	par(mfrow=c(1,2))

	if(exists("OUTDIR")){
		COND1<-paste0(cond1,".n",length(indS1))
		COND2<-paste0(cond2,".n",length(indS2))
		datout<-data.frame(gene=c(geneSymb2[indFC21],geneSymb2[indFC22]),
			status=rep(c(COND1,COND2),c(length(indFC21),length(indFC22))),
			COND1=COND1,
			COND2=COND2,
			len1=length(indS1),
			len2=length(indS2),
			numON1=c(flag1[indFC21],flag1[indFC22]),
			numON2=c(flag2[indFC21],flag2[indFC22]))

		prefix<-paste0("No.",formatC(compi,width=3,flag = "0"),".",COND1,"-vs-",COND2,".nDEG",nrow(datout))
		fnout=paste0(OUTDIR,"/FC.",prefix,".txt")
		cat("writing to: ",fnout,"\n")
		print(head(datout))
		write.table(datout,file=fnout,sep="\t",quote=F,row.names=F)
	}
}
VOLCANO<-function(vecIN,cond1="young",cond2="old"){
	if(compi==0){
		plot(0,0,type='n',axes=F,xlab="",ylab="")
		text(-1,-0.4,"There are three types of DEGs:
		(A) statistical DEGs: genes that have enough values in both groups
		(B) non-statistical DEGs by fold-change: genes where one group has enough data, but the other doesn't.
			But the other group still has at least one value. We detect the DEGs by fold-change cutoff.
		(C) non-statistical DEGs by on-off: genes where one group has no values at all, and the group has at least two values.

		-Dec 5, 2019         Peikai",pos =4,offset =0,xpd=T,cex=2/3)
		
	}
	compi<<-compi+1
	indS1<-which(vecIN==cond1)
	indS2<-which(vecIN==cond2)

	flag1<-rowSums(!is.na(data[,indS1]))
	flag2<-rowSums(!is.na(data[,indS2]))

	m1<-rowMeans(data[,indS1],na.rm=T)
	m2<-rowMeans(data[,indS2],na.rm=T)

	TH1<-sum(vecIN==cond1)/2
	TH2<-sum(vecIN==cond2)/2
	cat("TH1:",TH1,"TH2",TH2,"\n")
	par(mfrow=c(1,2))
	pMAT<-apply(data,1,function(x){
		y<-x[indS1]
		z<-x[indS2]
		y2<-y[!is.na(y)]
		z2<-z[!is.na(z)]
		if(length(y2)>TH1&length(z2)>TH2){
			res1<-t.test(y2,z2)
			return(c(res1$p.value,diff(res1$estimate),length(y2),length(z2)))
		}else{
			tmp1<-mean(z2)-mean(y2)
			return(c(2,tmp1,length(y2),length(z2)))
		}
	})
	print(table(pMAT[1,]<1,pMAT[2,]!=0))

	indT<-pMAT[1,]<1&pMAT[2,]!=0
	MAIN<-paste0("Comparison ",compi,"-A: ",cond1," (",sum(vecIN==cond1),"; blue)  vs. ",
		cond2," (",sum(vecIN==cond2),"; red)\n",
		"[",sum(indT)," genes were tested]")
	MAINb<-paste0("Comparison ",compi,"-B: ",cond1," (",sum(vecIN==cond1),"; blue)  vs. ",
		cond2," (",sum(vecIN==cond2),"; red)\n",
		"[",sum(indT)," genes were tested]")

	FDRq<-p.adjust(pMAT[1,indT])
	indSig<-which(FDRq<0.05)
	indSig2<-which(pMAT[1,indT]<0.05)
	#print(head(t(pMAT)))
	#plot(p1[2,indT],-log10(pMAT[1,indT]),pch=16,cex=2,col="grey",
	#main="all young vs all old")
	plot(pMAT[2,indT],-log10(FDRq),pch=".",cex=2,col="grey",
		xlab=paste0(cond2,"-",cond1),cex.main=2/3,
		ylab="-log10(FDR q)",
		main=MAIN)
	mtext(side=3,paste0("Blue: higher in ",cond1,";  Red: higher in ",cond2),cex=2/3)
	abline(v=0,lty=2)
	abline(h=0,lty=2)
	if(length(indSig)>0){
		points(pMAT[2,indT][indSig],-log10(FDRq[indSig]),pch=16,cex=2,
		col=hue_pal()(3)[-sign(pMAT[2,indT][indSig])+2])
		text(pMAT[2,indT][indSig],-log10(FDRq)[indSig],geneSymb2[indT][indSig],cex=2/3)
		DIFF<-pMAT[2,indT][indSig]
		
		#DRAWGENE(geneSymb2[indT][indSig][DIFF<0],geneSymb2[indT][indSig][DIFF>0])
		CAT_DEG2(geneSymb2[indT][indSig][DIFF>0],geneSymb2[indT][indSig][DIFF<0], BSIZE=5, WIDTH=60,Xpos=0.75)
	}

	par(mfrow=c(1,2))
	FDRq<-pMAT[1,indT]
	plot(pMAT[2,indT],-log10(pMAT[1,indT]),pch=".",cex=2,col="grey",
		xlab=paste0(cond2,"-",cond1),cex.main=2/3,
		ylab="-log10(p-values)",
		main=MAINb)
	mtext(side=3,paste0("Blue: higher in ",cond1,";  Red: higher in ",cond2),cex=2/3)
	abline(v=0,lty=2)
	abline(h=0,lty=2)

	if(length(indSig)>0){
		points(pMAT[2,indT][indSig],-log10(pMAT[1,indT][indSig]),pch=16,cex=2,
			col=hue_pal()(3)[-sign(pMAT[2,indT][indSig])+2])
		text(pMAT[2,indT][indSig],-log10(pMAT[1,indT])[indSig],geneSymb2[indT][indSig],cex=2/3)
		DIFF<-pMAT[2,indT][indSig]

#		DRAWGENE(geneSymb2[indT][indSig][DIFF<0],geneSymb2[indT][indSig][DIFF>0])
		CAT_DEG2(geneSymb2[indT][indSig][DIFF>0],geneSymb2[indT][indSig][DIFF<0], BSIZE=5, WIDTH=60,Xpos=0.75)
	}else{
		FDRq<-pMAT[1,indT]
		plot(pMAT[2,indT],-log10(pMAT[1,indT]),pch=".",cex=2,col="grey",
			xlab=paste0(cond2,"-",cond1),cex.main=2/3,
			ylab="-log10(p-values)",
			main=MAIN)
		mtext(side=3,paste0("Blue: higher in ",cond1,";  Red: higher in ",cond2),cex=2/3)
		abline(v=0,lty=2)
		abline(h=0,lty=2)

		text(pMAT[2,indT][indSig2],-log10(pMAT[1,indT])[indSig2],geneSymb2[indT][indSig2],cex=2/3)
	}
	if(exists("OUTDIR")&sum(indT[indSig])>0){
		COND1<-paste0(cond1,".n",length(indS1))
		COND2<-paste0(cond2,".n",length(indS2))
		datout<-data.frame(gene=geneSymb2[indT][indSig],
			status=c(cond1,cond2)[as.integer(pMAT[2,indT][indSig]>0)+1],
			COND1=COND1,
			COND2=COND2,
			len1=length(indS1),
			len2=length(indS2),
			numON1=flag1[indT][indSig],
			numON2=flag2[indT][indSig],
			t(pMAT[,indT][,indSig]))

		prefix<-paste0("No.",formatC(compi,width=3,flag = "0"),".",COND1,"-vs-",COND2,".nDEG",nrow(datout))
		fnout=paste0(OUTDIR,"/STATISTICAL.",prefix,".txt")
		cat("writing to: ",fnout,"\n")
		print(head(datout))
		write.table(datout,file=fnout,sep="\t",quote=F,row.names=F)
	}

##############
	FOLD_CHNAGE(vecIN,cond1,cond2,indT,indSig)
##############
	ONOFF(vecIN,cond1,cond2)

##############
	if(0){
	D1<-pMAT[4,]/2/TH2-pMAT[3,]/2/TH1
	D2<-pMAT[2,]
	par(mfrow=c(1,2))
	indO1<-D1>0.5&abs(D2)>1
	indO2<-D1< -0.5&abs(D2)>1
	oneside<-(pMAT[3,]>TH1|pMAT[4,]>TH2)
	plot(D1,D2,pch=16,col="lightgrey")
	points(D1[oneside],D2[oneside],pch=16,col="darkgrey")
	if(length(indSig)>0){
		points(D1[indT][indSig],D2[indT][indSig],pch=16,col=hue_pal()(3)[-sign(pMAT[2,indT][indSig])+2])
	}
	if(length(indSig2)>0){
		points(D1[indT][indSig2],D2[indT][indSig2],pch=16,col=hue_pal()(3)[-sign(pMAT[2,indT][indSig2])+2])
	}
	if(length(indO1)>0)text(D1[indO1],D2[indO1],geneSymb2[indO1],cex=1/2)
	if(length(indO2)>0)text(D1[indO2],D2[indO2],geneSymb2[indO2],cex=1/2)
	}

##############

##############
	res1<-list(pMAT=pMAT,indT=indT,indS1=indS1,indS2=indS2,indSig=indSig,indSig2=indSig2,cond1=cond1,cond2=cond2)
	return(res1)
}


