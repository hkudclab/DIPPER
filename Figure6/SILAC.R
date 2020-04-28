fn1<-list.files("march26",pattern="*txt",full.names=T)

LRES<-list()
for(i in 1:length(fn1)){
	tmpi<-read.csv(fn1[i],sep="\t")
	LRES[[i]]<-tmpi
}
names(LRES)<-gsub(".*IVD_|(_[cC]ombined.*txt)","",fn1)
####################################################################
load("../resources/MATRISOME.RData")
load("../resources/NABA.break.down.RData")
load("../resources/TF.CD.receptors.human.mouse.RData")

LNABA[["Collagens"]]<-as.character(MATRISOME[grep("^COL[0-9]",MATRISOME[,1]),1])
LNABA2<-LNABA
LNABA2[["TF"]]<-TF.human
LNABA2[["CD"]]<-surface.human
Lgenecat<-LNABA2[c(7,5,3,2,4,6)]


BARPLOT<-function(K=0,strin="Ratio"){
	NumHeavy<-sapply(LRES,function(X){
		indRatio<-grep("Ratio.*Heavy.*Light",colnames(X))
		sum(X[,indRatio+K]>5,na.rm=T)
	})
	Cat1<-gsub("OAF","AF",gsub("[0-9]+_.*","",names(NumHeavy)))
	Cat2<-factor(Cat1,levels=c( "NP_CDS", "AF_CDS", "NP_DS", "AF_DS"))
	par(mar=c(10,5,5,2));
	barplot(NumHeavy[order(Cat2)],main=strin,las=2)
}
MATPLOT<-function(K=0,strin="Ratio"){
	NumHeavy<-sapply(LRES,function(X){
		GN<-gsub(" .*","",gsub(".*GN=","",X[,5]))
		indRatio<-grep("Ratio.*Heavy.*Light",colnames(X))
		GNk<-GN[which(X[,indRatio+K]>5)]
		LENK<-sapply(Lgenecat,function(gCAT){
			length(intersect(gCAT,GNk))
		})
		LENK<-c(LENK,length(setdiff(GNk,unlist(Lgenecat))))
	})
	Cat1<-gsub("OAF","AF",gsub("[0-9]+_.*","",colnames(NumHeavy)))
	Cat2<-factor(Cat1,levels=c( "NP_CDS", "AF_CDS", "NP_DS", "AF_DS"))
	par(mar=c(10,5,5,2));
	barplot(NumHeavy[,order(Cat2)],col=c(seq(6),"grey"),main=strin,las=2)
	LEGEND<-rownames(NumHeavy);LEGEND[7]<-"Others"
	legend("topleft",
		fill=c(seq(6),"grey"),
		legend=LEGEND,cex=2/3)
	NumHeavy
}

pdf("SILAC.NumDetected.march27.pdf",width=12,height=4)
	par(mfrow=c(1,3))
	BARPLOT(1,"Light")
	BARPLOT(2,"Heavy")
	BARPLOT(0,"Ratio")

	MAT1<-MATPLOT(1,"Light")
	MAT2<-MATPLOT(2,"Heavy")
	MAT3<-MATPLOT(0,"Ratio")
dev.off()

pdf("SILAC.Ratio.Barcharts.march30.pdf",width=18)
	xMAX<-max(sapply(LRES,nrow))
	for(i in 1:8){
		par(mfrow=c(3,1))
		X<-LRES[[i]]
		GN<-gsub(" .*","",gsub(".*GN=","",X[,5]))
		KOLi<-sapply(GN,function(x){
			KOLx<-sapply(Lgenecat,function(gCAT){
				x%in%gCAT
			})
			which(KOLx)[1]
		})
		KOLi[is.na(KOLi)]<-"grey"

		indRatio<-grep("Ratio",colnames(X))

		bpp<-barplot(X[,indRatio],col=KOLi,border=0,xlim=c(0,xMAX),
			main=paste0("No. ",i," : ",names(LRES)[i],";      ",colnames(X)[indRatio]))
		text(bpp,-5,GN,srt=90,cex=1/3,pos=2,xpd=T,offset=0)

		bpp<-barplot(X[,indRatio+1],col=KOLi,border=0,xlim=c(0,xMAX),
			main=paste0("No. ",i," : ",names(LRES)[i],";      ",colnames(X)[indRatio+1]))
		text(bpp,-5,GN,srt=90,cex=1/3,pos=2,xpd=T,offset=0)

		bpp<-barplot(X[,indRatio+2],col=KOLi,border=0,xlim=c(0,xMAX),
			main=paste0("No. ",i," : ",names(LRES)[i],";      ",colnames(X)[indRatio+2]))
		text(bpp,-5,GN,srt=90,cex=1/3,pos=2,xpd=T,offset=0)
	}

dev.off()

#######################################################
LAF_CDS<-LRES[c(1,2,7,8)]
names(LAF_CDS)
LNP_CDS<-LRES[c(4,5)]
names(LNP_CDS)
LAF_DS<-LRES[3]
names(LAF_DS)
LNP_DS<-LRES[6]
names(LNP_DS)
CNT<-0
getAvrg<-function(LIN,K=0,fGAPDH=0){
	Lacc<-sapply(LIN,function(X){
		as.character(X[,"Accession"])
	})
	Lgn<-sapply(LIN,function(X){
		GN<-gsub(" .*","",gsub(".*GN=","",X[,5]))
	})
	UNQid<-unique(unlist(Lacc))
	UNQ_GN<-unlist(Lgn)[match(UNQid,unlist(Lacc))]

	KOLi<-sapply(UNQ_GN,function(x){
			KOLx<-sapply(Lgenecat,function(gCAT){
				x%in%gCAT
			})
			which(KOLx)[1]
		})
	KOLi[is.na(KOLi)]<-"grey"


	TMPi<-sapply(LIN,function(X){
		indRatio<-grep("Ratio.*Heavy.*Light",colnames(X))
		acci<-as.character(X[,"Accession"])
		X[match(UNQid,acci),indRatio+K]
		})
	rowM<-rowMeans(TMPi,na.rm=T)
	ORD<-order(rowM,decreasing=T)
	#print(rowM[ORD])
	indRatio<-grep("Ratio.*Heavy.*Light",colnames(LIN[[1]]))
	MAINstr<-paste0("No. ",CNT," : ",
		paste(names(LIN),collapse=";"),";      ",
		colnames(LIN[[1]])[indRatio+K],
		"    num=",length(rowM),
		"    numValid=",sum(!is.na(rowM)))
	MAINstr<-gsub("_7day","",MAINstr)
	MAINstr<-gsub("\\.\\.+","-",MAINstr)
	par(mar=c(2,4,2,2))
	bpp<-barplot(rowM[ORD],col=KOLi[ORD],border=0,xlim=c(0,500),
		main=MAINstr)
	text(bpp,5,UNQid[ORD],srt=90,cex=1/3,pos=4,xpd=T,offset=0)
	text(bpp,-5,UNQ_GN[ORD],srt=90,cex=1/3,pos=2,xpd=T,offset=0)
	abline(v=bpp[seq(1,500,by=100)],lty=3,col=3,lwd=2)
	box("plot")
	if(fGAPDH==1){
		indGAPDH<-which(UNQ_GN[ORD]=="GAPDH")
		strN<-paste0(UNQ_GN[ORD][seq(indGAPDH)],"=",UNQid[ORD][seq(indGAPDH)])
		print(length(ORD))
	}else if(fGAPDH==2){
		indGAPDH<-length(which(rowM>1))
		print(indGAPDH)
		strN<-paste0(UNQ_GN[ORD][seq(indGAPDH)],"=",UNQid[ORD][seq(indGAPDH)])
		print(length(ORD))
	}else if(fGAPDH==3){
		indGAPDH<-length(which(rowM>1))
		strN<-paste0(UNQ_GN[ORD][seq(indGAPDH)],"=",UNQid[ORD][seq(indGAPDH)])
		strN<-data.frame(strN,level=rowM[ORD][seq(indGAPDH)])
	}
	return(strN)
}



pdf("SILAC.averge.Ratio.Barcharts.march30.pdf",width=18)
	par(mfrow=c(4,1))
	Rt_CDS_AF<-getAvrg(LAF_CDS,0)
	Rt_CDS_AF<-getAvrg(LAF_DS,0)
	Rt_DS_NP<-getAvrg(LNP_CDS,0)
	Rt_DS_NP<-getAvrg(LNP_DS,0)

	par(mfrow=c(4,1))
	Lt_CDS_AF<-getAvrg(LAF_CDS,1)
	Lt_DS_AF<-getAvrg(LAF_DS,1)
	Lt_CDS_NP<-getAvrg(LNP_CDS,1)
	Lt_DS_NP<-getAvrg(LNP_DS,1)

	par(mfrow=c(4,1))
	Hv_CDS_AF<-getAvrg(LAF_CDS,2)
	Hv_DS_AF<-getAvrg(LAF_DS,2)
	Hv_CDS_NP<-getAvrg(LNP_CDS,2)
	Hv_DS_NP<-getAvrg(LNP_DS,2)

dev.off()

library(gplots)
L4_Lt<-list(CDS_AF=Lt_CDS_AF,
	CDS_NP=Lt_CDS_NP,
	DS_AF=Lt_DS_AF,
	DS_NP=Lt_DS_NP)

L4_Hv<-list(CDS_AF=Hv_CDS_AF,
	CDS_NP=Hv_CDS_NP,
	DS_AF=Hv_DS_AF,
	DS_NP=Hv_DS_NP)


	Lt_CDS_AF<-getAvrg(LAF_CDS,1,fGAPDH=3)
	Lt_DS_AF<-getAvrg(LAF_DS,1,fGAPDH=3)
	Lt_CDS_NP<-getAvrg(LNP_CDS,1,fGAPDH=3)
	Lt_DS_NP<-getAvrg(LNP_DS,1,fGAPDH=3)
L4_Lt<-list(CDS_AF=Lt_CDS_AF,
	CDS_NP=Lt_CDS_NP,
	DS_AF=Lt_DS_AF,
	DS_NP=Lt_DS_NP)

	Hv_CDS_AF<-getAvrg(LAF_CDS,2,fGAPDH=3)
	Hv_DS_AF<-getAvrg(LAF_DS,2,fGAPDH=3)
	Hv_CDS_NP<-getAvrg(LNP_CDS,2,fGAPDH=3)
	Hv_DS_NP<-getAvrg(LNP_DS,2,fGAPDH=3)
L4_Hv<-list(CDS_AF=Hv_CDS_AF,
	CDS_NP=Hv_CDS_NP,
	DS_AF=Hv_DS_AF,
	DS_NP=Hv_DS_NP)
save(L4_Lt,L4_Hv,Lgenecat,file="SILAC.heavy.light.pooled.RData")

source("step5.help.funcs.R")

pdf("Light.Heavy.above.GAPDH.venn.apr-24.pdf")
	v1<-venn(L4_Lt)
	inter<-attr(v1,"intersections")
	inter2<-sapply(inter,function(x)gsub("=.*","",x))
	showVenn(inter,BSIZE=5,WD=130)
	showVenn(inter2,BSIZE=5,WD=130)

###############
	v1<-venn(L4_Hv)
	inter<-attr(v1,"intersections")
	inter2<-sapply(inter,function(x)gsub("=.*","",x))
	showVenn(inter,BSIZE=5,WD=130)
	showVenn(inter2,BSIZE=5,WD=130)
dev.off()
