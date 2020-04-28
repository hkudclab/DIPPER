
oydata<-read.csv("old-and-young.aug21-2019.txt",sep="\t")

head(oydata)
geneSymb<-as.character(oydata[,3])
geneSymb2<-gsub(";.*$","",geneSymb)

setdiff(geneSymb2,annot[,"gene_name"])
nonMatrisome<-setdiff(geneSymb2,MATRISOME[,1])

data0<-as.matrix(oydata[,-seq(3)])

agegrp<-"old_young"
sampleIDS<-colnames(data0)
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
load("../resources/MATRISOME.RData")
load("../resources/NABA.break.down.RData")

LNABA[["Collagens"]]<-as.character(MATRISOME[grep("^COL[0-9]",MATRISOME[,1]),1])
LNABA2<-LNABA
Lgenecat<-LNABA2[c(7,5,3,2,4,6)]
Lgenecat[["nonMatrisome"]]<-setdiff(geneSymb2,unlist(Lgenecat))


save(oydata,geneSymb2,
	data,
	sampleIDS,
	Ages,Compartments,Directions,Levels,
	Lgenecat,
	file="proteome.young-old.n66.RData")


