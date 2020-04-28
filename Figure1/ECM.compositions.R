library(scales)
library(stringr)

load("../preprocess/proteome.young-old.n66.RData")

f1<-factor(Compartments,levels=c("NP", "IAF_NP", "IAF", "OAF"))
ORD1<-order(paste0(Ages,".",as.integer(f1),".",Directions,".",Levels))

for(i in 1:length(Lgenecat)){
	LgeneMATxy<-apply(data[,ORD1],2,function(y){
		intersect(geneSymb2[which(y>0)],Lgenecat[[i]])
		})
	cat(names(Lgenecat)[i],":",length(unique(unlist(LgeneMATxy))),"\n")
	#print(sort(table(unlist(LgeneMATxy)),decreasing=T))
}


MATxy<-t(apply(data[,ORD1],2,function(x){
	sapply(Lgenecat,function(y){
		length(intersect(geneSymb2[which(x>0)],y))
		})
	}))
mat1<-t(t(sapply(Lgenecat,function(y){
		length(intersect(geneSymb2[rowSums(!is.na(data))==66],y))
		})))
t(t(sapply(Lgenecat,function(y){
		paste0(sort(intersect(geneSymb2[rowSums(!is.na(data))==66],y)),collapse=", ")
		})))
RATIOxy<-apply(MATxy,1,function(x)x/sum(x))

t(t(sapply(Lgenecat,function(y){
		paste0(sort(intersect(geneSymb2[rowSums(!is.na(data))==66],y)),collapse=", ")
		})))
    
##############################################
library(grid)
library(lattice)
library(gridExtra)
library(stringr)
library(gplots)

ADDTITLE<-function(str1="",Xpos=0.76,BSIZE=5,WD=90){
	layout(t(matrix(seq(2))),widths=c(6,4))
	v1<-venn(Lij)
	mytheme <- gridExtra::ttheme_default(base_size=BSIZE,
		padding = unit(c(1, 3), "mm"))

	CNT<<-CNT+1

	USR<-par('usr')

	x<-mean(USR[1:2])
	y<-diff(USR[3:4])*1+USR[3]
	text(x,y,paste0("NO.",CNT,": ",str1,"\ntotal:",length(unique(unlist(Lij)))),xpd=T)

	INTER<-attr(v1,"intersections")
	MAT<-as.matrix(sapply(INTER[sapply(INTER,length)<80],function(x){
		str_wrap(paste0(paste0(x,collapse=", "),"\n(n=",length(x),")"),WD)
		}))
	if(nrow(MAT)>0){
		pushViewport(viewport(x=Xpos,y=0.5))
		grid.table(MAT,theme=mytheme )
		popViewport()
	}
}
##############################################
COMPART<-unique(Compartments)

pdf("ECM.compositions.venn.Mar20.pdf",height=6,width=10)
	CNT<-0
	Lij<-sapply(COMPART,function(y){
		indy<-which(Compartments==y)
		geneSymb2[rowSums(!is.na(data[,indy]))>0]
		})
	ADDTITLE("young & aged, all genes")

	for(j in 1:length(Lgenecat)){
		Lij<-sapply(COMPART,function(y){
			indy<-which(Compartments==y)
			intersect(geneSymb2[rowSums(!is.na(data[,indy]))>0],Lgenecat[[j]])
			})
		ADDTITLE(paste0("young&aged, ",names(Lgenecat)[j]))
	}

	Lij<-sapply(COMPART,function(y){
		indy<-which(Compartments==y&Ages=="young")
		geneSymb2[rowSums(!is.na(data[,indy]))>0]
		})
	ADDTITLE("young only, all genes")

	Lij<-sapply(COMPART,function(y){
		indy<-which(Compartments==y&Ages=="old")
		geneSymb2[rowSums(!is.na(data[,indy]))>0]
		})
	ADDTITLE("aged only, all genes")

	for(i in 1:2){
		myAge<-c("young","old")[i]

		for(j in 1:length(Lgenecat)){
			Lij<-sapply(COMPART,function(y){
				indy<-which(Compartments==y&Ages==myAge)
				intersect(geneSymb2[rowSums(!is.na(data[,indy]))>0],Lgenecat[[j]])
				})
			ADDTITLE(paste0(myAge," only, ",names(Lgenecat)[j]))

		}
	}
dev.off()


