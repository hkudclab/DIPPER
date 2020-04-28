library(scales)

load("../resources/anita.RData")
load("../resources/MATRISOME.RData")

oydata<-read.csv("old-and-young.aug21-2019.txt",sep="\t")

head(oydata)
geneSymb<-as.character(oydata[,3])
geneSymb2<-gsub(";.*$","",geneSymb)

genesAnita<-gsub(" .*","",annots[,2])
indAnita1<-which(genesAnita%in%MATRISOME[,1])
indAnita2<-which(genesAnita%in%CORE[,1])

indAnita3<-match(intersect(geneSymb2,genesAnita),genesAnita)
indAnita4<-match(intersect(geneSymb2,MATRISOME[,1]),genesAnita)
library(scales)

boxplot(dat.log,las=2)
#pdf("anita.scatter.plot.pdf",width=9,height=9)
#	plot(data.frame(dat.log),pch=".")
#dev.off()

PID<-substr(colnames(dat.log),1,3)
GROUP<-substr(PID,1,1)
AFNP<-substr(colnames(dat.log),4,5)

AGEComp<-paste0(GROUP,"_",AFNP)
DATmean<-sapply(unique(AGEComp),function(x){
	rowMeans(dat.log[,AGEComp%in%x])
})

KOLanita<-rev(hue_pal()(2))[as.integer(as.factor(GROUP))]
PCHanita<-c(18,16)[as.integer(as.factor(AFNP))]

res1<-prcomp(t(dat.log))
pc1<-paste0(paste0("PC",seq(8)," ("),percent(res1$sdev^2/sum(res1$sdev^2)),")")

source("step9.help.funcs.SVM.R")
SVMPLOT(res1,GROUP,pc1,DIM="12")
SVMPLOT(res1,AFNP,pc1,DIM="12")

plot(res1$x[,1],res1$x[,2],
	pch=PCHanita,
	cex=4,xlab=pc1[1],ylab=pc1[2],
	col=KOLanita)
legend("bottom",fill=hue_pal()(4),
	ncol=2,
	legend=levels(as.factor(PID)),cex=1.5)
text(res1$x[,1],res1$x[,2],colnames(dat.log))

##################
young<-dat.log[,GROUP=="C"]
old<-dat.log[,GROUP=="D"]
yAF<-dat.log[,c(1,3)]
yNP<-dat.log[,c(2,4)]
oAF<-dat.log[,c(5,7)]
oNP<-dat.log[,c(6,8)]
yAFNPdiff<-rowMeans(yAF)-rowMeans(yNP)
oAFNPdiff<-rowMeans(oAF)-rowMeans(oNP)
oyAFdiff<-rowMeans(oAF)-rowMeans(yAF)
oyNPdiff<-rowMeans(oNP)-rowMeans(yNP)
##################
plot(as.dendrogram(hclust(as.dist(1-cor(dat.log,method="pearson")))))

plot(as.dendrogram(hclust(as.dist(1-cor(dat.log[indAnita1,],method="pearson")))))
plot(as.dendrogram(hclust(as.dist(1-cor(dat.log[indAnita2,],method="pearson")))))
plot(as.dendrogram(hclust(as.dist(1-cor(dat.log[indAnita3,],method="pearson")))))
plot(as.dendrogram(hclust(as.dist(1-cor(dat.log[indAnita4[!is.na(indAnita4)],],method="spearman")))))


MAPLOT<-function(dat1,dat2,str1,str2){
	x<-rowMeans(dat1)
	y<-rowMeans(dat2)
	ind1<-which(abs(x-y)>3&(x+y)>10&genesAnita!="")
	print(genesAnita[ind1])
	COL<-c("#F8766D", "#619CFF")[as.integer((x-y)[ind1]>0)+1]
	Leg1<-paste0(str1,":",sum(COL=="#619CFF"))
	Leg2<-paste0(str2,":",sum(COL=="#F8766D"))

	plot(x,y,pch=".",xlab=str1,ylab=str2)
	points(x[ind1],y[ind1],col=COL,cex=2,pch=16)
	text(x[ind1],y[ind1],genesAnita[ind1])
	legend("topleft",legend=c(Leg2,Leg1),fill=c("#F8766D", "#619CFF"))

	CNT<<-CNT+1
	outDAT<-data.frame(genes=genesAnita,x,y)[ind1,]
	colnames(outDAT)[2:3]<-c(str1,str2)
	#fnout<-paste0("Anita/NO.",CNT,".",str2,"-vs-",str1,".txt")
	#write.table(outDAT,file=fnout,quote=F,row.names=F,sep="\t")
	outDAT
}
pdf("Anita/anita.MAPLOT.apr22.pdf",width=9,height=9)
	CNT<-0
	res1<-MAPLOT(yAF,yNP,"yAF","yNP")
	res2<-MAPLOT(oAF,oNP,"oAF","oNP")
	res3<-MAPLOT(yAF,oAF,"yAF","oAF")
	res4<-MAPLOT(yNP,oNP,"yNP","oNP")

	MAPLOT(young,old,"young","old")

dev.off()
##################################
library(gplots)
load("../SILAC/SILAC.YND.AGD.ECM.RData")
AGDAF<-names(which(DATAF[2,]>0))
AGDNP<-names(which(DATNP[2,]>0))

L4<-split(as.character(res4[,1]),factor(res4$oNP>res4$yNP))
L4[["AGDNP"]]<-AGDNP
venn(L4)

L3<-split(as.character(res3[,1]),factor(res3$oAF>res3$yAF))
L3[["AGDAF"]]<-AGDAF
venn(L3)

ind3<-which(genesAnita%in%AGDAF)
boxplot(dat.log[ind3,],las=2)
boxplot(cbind(rowMeans(yAF),rowMeans(oAF))[ind3,],las=2)

ind4<-which(genesAnita%in%AGDNP)
boxplot(cbind(rowMeans(yNP),rowMeans(oNP))[ind4,],las=2)

##################################
library(gplots)

fnDEG<-list.files(path="Anita",pattern="NO.[3-5].*txt",full.names=T)
LDEG<-rep(list(),3)
for(i in 1:3){
	tmpi<-read.csv(fnDEG[i],sep="\t")
	LDEG[[i]]<-tmpi
}
names(LDEG)<-c("AF","NP","both")
Lyoung<-sapply(LDEG,function(X){
		unique(as.character(X[X[,2]>X[,3],1]))
	})
Lold<-sapply(LDEG,function(X){
		unique(as.character(X[X[,2]<X[,3],1]))
	})

pdf("Anita/Venn.pdf")
	v1<-venn(Lyoung)
	attr(v1,"intersections")
	v2<-venn(Lold)
	attr(v2,"intersections")

dev.off()

