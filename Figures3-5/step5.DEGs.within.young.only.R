

library(scales)
library(stringr)

load("../resources/MATRISOME.RData")
load("../resources/NABA.break.down.RData")
load("../resources/gencode.v25.annotation.gene.gtf.RData")
load("../resources/TF.CD.receptors.human.mouse.RData")


oydata<-read.csv("old-and-young.aug21-2019.txt",sep="\t")

load("../preprocess/proteome.young-old.n66.RData")


data0<-data
Lpheno<-strsplit(sampleIDS,"\\.")
Ages<-tolower(sapply(Lpheno,function(x)x[2]))

data<-data0[,Ages%in%unlist(strsplit(agegrp,"_"))]
sampleIDS<-colnames(data)
Lpheno<-strsplit(sampleIDS,"\\.")

Levels<-sapply(Lpheno,function(x)x[1])
Ages<-tolower(sapply(Lpheno,function(x)x[2]))
Directions<-sapply(Lpheno,function(x)x[3])
Compartments<-sapply(Lpheno,function(x)x[4])
Compartments[is.na(Compartments)]<-Directions[is.na(Compartments)]
Compartments2<-Compartments
Compartments2[grep("AF",Compartments)]<-"allelse"
Compartments2[Compartments=="IAF_NP"]<-"NPoriN"
Compartments2[Compartments=="NP"]<-"NPorNi"

Compartments3<-Compartments
Compartments3[grep("AF",Compartments)]<-"others"

AP<-Directions
AP[!AP%in%c("A","P")]<-NA
LR<-Directions
LR[!LR%in%c("L","R")]<-NA
#######################################
source("step5.helper.funcs.R")
Lgenecat<-LNABA2[c(7,5,3,2,4,6)]
Lgenecat[["nonMatrisome"]]<-setdiff(geneSymb2,unlist(Lgenecat))

#######################################
pdf("DEG.young.only.DEC4.pdf",width=12)
	compi<-0

	isOAF<-Compartments
	isOAF[isOAF!="OAF"]<-"NONoAF"
	vec1<-paste0(Ages,"_",isOAF)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_NONoAF","young_OAF")
	withANITA("young AF - young NP (microarray)",p10)
#dev.off()
#exit
##############
	isOAF<-Compartments
	isOAF[isOAF!="OAF"&isOAF!="IAF"]<-"NONoiAF"
	isOAF[isOAF=="OAF"|isOAF=="IAF"]<-"oiAF"
	vec1<-paste0(Ages,"_",isOAF)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_NONoiAF","young_oiAF")
	withANITA("young AF - young NP (microarray)",p10)
##############
	isOAF<-Compartments
	isOAF[isOAF!="IAF"]<-"NONiAF"
	vec1<-paste0(Ages,"_",isOAF)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_NONiAF","young_IAF")
##############
	isLower<-Levels
	isLower[isLower!="L5_S1"]<-"upper"
	vec1<-paste0(Ages,"_",Levels)
	vec2<-paste0(Ages,"_",isLower)
	table(vec1,Ages)
	table(vec2,Ages)

	p10<-VOLCANO(vec1,"young_L5_S1","young_L4_5")
	p10<-VOLCANO(vec1,"young_L5_S1","young_L3_4")
	p10<-VOLCANO(vec1,"young_L3_4","young_L4_5")
	p10<-VOLCANO(vec2,"young_L5_S1","young_upper")
##############
	LRAP2<-rep(NA,length(LR))
	LRAP2[!is.na(LR)]<-"LR"
	LRAP2[!is.na(AP)]<-"AP"
	vec1<-paste0(Ages,"_",Compartments,"_",LRAP2)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_OAF_AP","young_OAF_LR")
	p10<-VOLCANO(vec1,"young_IAF_NP_AP","young_IAF_NP_LR")

##############
	ioAF<-Compartments
	ioAF[(ioAF!="OAF"&ioAF!="IAF")|is.na(LR)]<-"others"
	vec1<-paste0(Ages,"_",ioAF,"_LR")
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_LR","young_OAF_LR")
##############
	iNP<-Compartments
	iNP[iNP!="OAF"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_NP","young_OAF")
##############
	iNP<-Compartments
	iNP[iNP!="OAF"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP,"_",c("nonLR","LR")[as.integer(!is.na(LR))+1])
	table(vec1,LR)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_NP_LR","young_OAF_LR")
##############
	iNP<-Compartments
	iNP[iNP!="OAF"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP,"_",c("nonAP","AP")[as.integer(!is.na(AP))+1])
	table(vec1,AP)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_IAF_NP_AP","young_OAF_AP")
##############
	iNP<-Compartments
	iNP[iNP!="NP"&iNP!="IAF_NP"]<-"others"
	vec1<-paste0(Ages,"_",iNP)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_NP","young_IAF_NP")
##############
	vec1<-paste0(Ages,"_",Compartments)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_IAF","young_OAF")

##############
	vec1<-paste0(Ages,"_",LR)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_L","young_R")
##############
	vec1<-paste0(Ages,"_",LR,"_",Compartments!="IAF")
	table(vec1,Compartments)
	p10<-VOLCANO(vec1,"young_L_TRUE","young_R_TRUE")
##############
	LR2<-LR;LR2[!is.na(LR2)]<-"LR"
	AP2<-AP;AP2[!is.na(AP2)]<-"AP"
	vec1<-paste0(Ages,"_",LR2,AP2)
	table(vec1,Ages)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_LRNA","young_NAAP")
############## All A vs all P
	vec1<-paste0(Ages,"_",AP)
	table(vec1,Ages)
	p10<-VOLCANO(vec1,"young_A","young_P")
##############
	LR2<-LR;LR2[!is.na(LR2)&Compartments!="IAF"]<-"LR"
	AP2<-AP;AP2[!is.na(AP2)&Compartments!="IAF"]<-"AP"
	vec1<-paste0(Ages,"_",LR2,AP2)
	table(vec1,Compartments)
	split(sampleIDS,vec1)
	p10<-VOLCANO(vec1,"young_LRNA","young_NAAP")
dev.off()
