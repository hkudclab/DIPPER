load("MRI.data.for.old.apr26.RData")
####################################################################
load("../resources/MATRISOME.RData")
load("../resources/NABA.break.down.RData")

LNABA[["Collagens"]]<-as.character(MATRISOME[grep("^COL[0-9]",MATRISOME[,1]),1])
LNABA2<-LNABA
Lgenecat<-LNABA2[c(7,5,3,2,4,6)]
####################################################################
oydata<-read.csv("../old-and-young/old-and-young.aug21-2019.txt",sep="\t")
data0<-as.matrix(oydata[,-seq(3)])

head(oydata)
geneSymb<-as.character(oydata[,3])
geneSymb2<-gsub(";.*$","",geneSymb)

agegrp<-"old_young"
sampleIDS<-colnames(data0)
Lpheno<-strsplit(sampleIDS,"\\.")
Ages<-tolower(sapply(Lpheno,function(x)x[2]))


data<-data0[,Ages%in%unlist(strsplit(agegrp,"_"))]
sampleIDS<-colnames(data)
Lpheno<-strsplit(sampleIDS,"\\.")

Levels<-sapply(Lpheno,function(x)x[1])
Ages<-tolower(sapply(Lpheno,function(x)x[2]))
APLRC<-sapply(Lpheno,function(x)x[3])
AFNP<-sapply(Lpheno,function(x)x[4])
AFNP[is.na(AFNP)]<-APLRC[is.na(AFNP)]
AFNP2<-AFNP
AFNP2[grep("AF",AFNP)]<-"allelse"
AFNP2[AFNP=="IAF_NP"]<-"NPoriN"
AFNP2[AFNP=="NP"]<-"NPorNi"

AFNP3<-AFNP
AFNP3[grep("AF",AFNP)]<-"others"

AP<-APLRC
AP[!AP%in%c("A","P")]<-NA
LR<-APLRC
LR[!LR%in%c("L","R")]<-NA

AGE_AFNP<-paste0(Ages,"_",AFNP)
#######################################
ind12<-match(MRIL,colnames(data))
KOLI<-sapply(geneSymb2,function(x){
	vec<-sapply(Lgenecat,function(cat){
		x%in%cat
	})
	which(vec)[1]
})
indECM<-which(KOLI%in%seq(6))
corAll<-cor(MRIm,t(data[,ind12]),use="pairwise.complete.obs")[1,]
corAllpval<-sapply(seq(3100),function(i){
	val<-tryCatch({
		resi<-cor.test(MRIm,data[i,ind12],use="pairwise.complete.obs")
		resi$p.value
	},  error = function(e){
		NA
	})	
	val
})
table(!is.na(corAllpval))
indSignif<-which(corAllpval<0.05)
str(indSignif)
indSignifECM<-intersect(indECM,indSignif)

#indSignifECM<-indSignif
indPos<-intersect(which(corAll>0),indSignifECM)
indNeg<-intersect(which(corAll<0),indSignifECM)

ORD1<-indPos[order(corAll[indPos])]
ORD1<-ORD1[order(KOLI[ORD1])]

ORD2<-indNeg[order(corAll[indNeg])]
ORD2<-ORD2[order(KOLI[ORD2])]

ORD<-c(ORD2,ORD1)

pdf("ECM.corr.with.MRI.all.three.levels.pdf",width=12)
	bpp<-barplot(corAll[ORD],col=KOLI[ORD])
	abline(h=c(-0.5,0.5),lty=2,col=3)
	text(bpp,sign(corAll[ORD])*-0.02,geneSymb2[ORD],
		srt=90,pos=c(4,2)[as.integer(corAll[ORD]>0)+1],
		offset=0)

dev.off()
##################################################
library("glmnet")
library("mvtnorm") 

load("MRI.data.for.young.apr27.RData")
ind13<-match(gsub("old","YOUNG",MRIL),toupper(colnames(data)))
table(toupper(yMRIL)==toupper(colnames(data))[ind13])

KeyECMs<-geneSymb2[ORD]
library(mice)

trainMat<-t(data[ORD,ind12])
colnames(trainMat)<-KeyECMs
testMat<-t(data[ORD,ind13])
colnames(testMat)<-KeyECMs


#tempData <- mice(trainMat,m=5,maxit=50,meth='pmm',seed=500)
#tempData2 <- mice(testMat,m=5,maxit=50,meth='pmm',seed=500)
#save(tempData ,tempData2, file="Imputed.for.MRI.prediction.Lasso.RData")
#load("Imputed.for.MRI.prediction.Lasso.RData")
completedDataTrain <- complete(tempData,1)
completedDataTest <- complete(tempData2 ,1)
indSEL<-which(colSums(is.na(completedDataTrain))==0&colSums(is.na(completedDataTest ))==0)
completedDataTrain2<-completedDataTrain[,indSEL]
colnames(completedDataTrain2)<-KeyECMs[indSEL]
completedDataTest2<-completedDataTest[,indSEL]
colnames(completedDataTest2)<-KeyECMs[indSEL]

#trainDat[is.na(trainDat)]<-0


#######################
trainDat<-data.frame(MRIm,completedDataTrain2)
fit1<-lm(MRIm~.,data=trainDat)
summary(fit1)

mylogit <- glm(MRIm~ ., data = trainDat, family = "gaussian")
summary(mylogit )

#######################
lambda <- 10^seq(10, -2, length = 100)
Xtrain<-as.matrix(completedDataTrain2)
lasso.mod <- glmnet(Xtrain, MRIm , alpha = 1, lambda = lambda)
cv.out <- cv.glmnet(Xtrain, MRIm , alpha = 1)
bestlam<-cv.out$lambda.min



MRIpred<-predict(fit1,completedDataTest2)
plot(yMRIm,MRIpred)

lasso.pred <- predict(lasso.mod, s = bestlam, newx = as.matrix(completedDataTest2))
mean((lasso.pred-yMRIm)^2)

abline(0,1,col=2)
CompartTest<-factor(AFNP[ind13],levels=rev(unique(AFNP)))
boxplot(yMRIm~CompartTest,xlab="compartments",ylab="young MRI intensities")
boxplot(lasso.pred~CompartTest,xlab="compartments",ylab="LASSO predicted")

library(ggplot2)

dat34<-data.frame(yMRIm=yMRIm,
	"lasso"=lasso.pred,
	compartments=CompartTest)
colnames(dat34)[2]<-"lasso"

library(scales)
COL_AFNP<-hue_pal()(4)[as.integer(factor(AFNP))]
PCH_ages<-c(17,16)[as.integer(factor(Ages))]
COL_levels<-seq(3)[as.integer(factor(Levels))]

##########################
pdf("Lasso.prediction.of.young.MRI.pdf")
	plot(yMRIm,lasso.pred,
		xlab="real young MRI intensities",
		cex=3,pch=16,col=COL_AFNP[ind13],
		ylab="LASSO predicted young MRI intensities")
	mtext(paste0("r=",signif(cor(yMRIm,lasso.pred))),cex=2)
	text(yMRIm,lasso.pred,gsub("[yY]oung.","",sampleIDS[ind13]),xpd=T,cex=2/3)
	legend("bottomright",fill=hue_pal()(4),
		legend=levels(factor(AFNP)))

	p34<- ggplot(dat34, aes(x=CompartTest, y=yMRIm)) + 
		geom_violin( scale = "width")+ 
		geom_boxplot(width=0.4,outlier.shape =NA)+ 
		geom_jitter(shape=16, size=4,color="#606060",position=position_jitter(0.2)) +
		ggtitle("real MRI")
	plot(p34)
	p34<- ggplot(dat34, aes(x=CompartTest, y=lasso)) + 
		geom_violin( scale = "width")+ 
		geom_boxplot(width=0.4,outlier.shape =NA)+ 
		geom_jitter(shape=16, size=4,color="#606060",position=position_jitter(0.2)) +
		ggtitle("lasso predicted MRI")
	plot(p34)

######################################
#         LASSO
######################################
  library(pROC)

  df<-data.frame(predictions=as.numeric(lasso.pred),labels=grepl("OAF",yMRIL))

  pROC_obj <- roc(df$labels,df$predictions,
            smoothed = TRUE,
            ci=TRUE, ci.alpha=0.9, stratified=FALSE,
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)


  pROC_obj <- roc(df$labels,df$predictions,
            smoothed = TRUE,
            ci=TRUE, ci.alpha=0.9, stratified=FALSE,
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)
  sens.ci <- ci.se(pROC_obj)
  plot(sens.ci, type="shape", col="lightblue")
## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.
  plot(sens.ci, type="bars")

dev.off()




