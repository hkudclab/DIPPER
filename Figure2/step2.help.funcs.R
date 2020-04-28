#LNABA[["CORE"]]<-as.character(CORE[,1])
#LNABA[["MATRISOME"]]<-as.character(MATRISOME[,1])
LNABA[["Collagens"]]<-as.character(MATRISOME[grep("^COL[0-9]",MATRISOME[,1]),1])

LNABA2<-LNABA
LNABA2[["TF"]]<-TF.human
LNABA2[["CD"]]<-surface.human


library(grid)
library(lattice)
library(gridExtra)
library(stringr)
##########################################
library(e1071)
SVMPLOT<-function(PCAi,DIM="12",indS,i){
	PCH_ages2<-PCH_ages[indS]
	COL_AFNP2<-COL_Compartments[indS]
	COL_levels2<-COL_levels[indS]
	COL_IO2<-COL_IO[indS]
	percenti<-Lpercent[[i]]

	KOLi<-cbind(COL_AFNP2,COL_IO2,COL_levels2)[,i]
	Y<-data.frame(Ages,innerOuter,gsub("L3_4|L4_5","Upper",Levels))[indS,i]
	TEXT<-c("Ages","inner/outer","Levels")[i]
	if(length(unique(Y))<2)return()
	if(DIM=="12"){
		inddim<-c(1,2)
	}else{
		inddim<-c(1,3)
	}
	print(length(unique(Y)))
	print(table(Y))
	print(dim(PCAi))
	dat = data.frame(PCAi$x[,inddim], y = Y)
	svmfit = svm(y ~ ., data = dat, kernel = "polynomial",coef0	=1)
	#plot(svmfit, dat)
	px1<-seq(range(PCAi$x[,inddim[1]])[1]*1.05,range(PCAi$x[,inddim[1]])[2]*1.05,len=20)
	px2<-seq(range(PCAi$x[,inddim[2]])[1]*1.05,range(PCAi$x[,inddim[2]])[2]*1.05,len=20)

	if(DIM=="12"){
		xgrid = expand.grid(PC1 = px1, PC2= px2)
	}else{
		xgrid = expand.grid(PC1 = px1, PC3= px2)
	}

	if(1){

		ygrid = predict(svmfit, xgrid)

		func = predict(svmfit, xgrid, decision.values = TRUE)
		func = attributes(func)$decision

		plot(PCAi$x[,inddim],cex=2,xlab=percenti[1],ylab=percenti[2],
			main=TEXT,
			pch=PCH_ages2,
			col=KOLi)
		contour(px1, px2, matrix(func, 20, 20), level = 0, add = TRUE)
		contour(px1, px2, matrix(func, 20, 20), level = 0.5, add = TRUE, col = "red", lwd = 2,lty=2)
		contour(px1, px2, matrix(func, 20, 20), level = -0.5, add = TRUE, col = "blue", lwd = 2,lty=2)

	}
}

PCAfunc<-function(indG1,IMPUTE,indi,NUMsample,indS=seq(66),SVMage=F,SVMlevel=F,SVMio=F){
	COL_AFNP2<-COL_Compartments[indS]
	PCH_ages2<-PCH_ages[indS]
	COL_levels2<-COL_levels[indS]

	i<-2*indi-1
	j<-2*indi
	PCAi<-LPCA[[i]]
	percenti<-Lpercent[[i]]
	PCAj<-LPCA[[j]]
	percentj<-Lpercent[[j]]


	geneSymb3<-geneSymb2[indG1]
	MAINSTR<-paste0("PCA on samples ",IMPUTE,"\nBased on ",length(geneSymb3),
		" genes with >",NUMsample," valid values")

	plot(1,type='n',axes=F,xlab="",ylab="")
	text(1,1,
		paste0("These are the ",length(geneSymb3),
		" genes that are detected in >",NUMsample," samples:\n\n\n",
		str_wrap(paste0(sort(geneSymb3),collapse=", ")),200),xpd=T,cex=1/2)

	mytheme <- gridExtra::ttheme_default(base_size =8,padding = unit(c(1, 1), "mm"))

	L2<-sapply(LNABA2,function(x)intersect(x,geneSymb3))
	L2[["others"]]<-setdiff(geneSymb3,unlist(L2))
	names(L2)<-paste0(c(names(LNABA2),"others")," (",sapply(L2,length),")")
	str2<-sapply(L2,paste0,collapse=", ")
	mat1<-as.matrix( str_wrap(str2,60))

	plot.new()
	mytheme <- gridExtra::ttheme_default(base_size =6,padding = unit(c(1, 1), "mm"))
	rownames(mat1)<-names(L2)
	grid.table(mat1,theme=mytheme )


	pci<- as.integer(gsub(".*\\(|%.*$","",percenti))
	bpp<-barplot(pci,main="Top PCs",ylab="variance captured per PC (%)")
	text(bpp,0,percenti,srt=90,offset=0,cex=1/2,pos=2,xpd=T)


	plot(PCAi$x[,1:2],cex=2,xlab=percenti[1],ylab=percenti[2],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_AFNP2)
	legend("topleft",fill=hue_pal()(4),
		legend=levels(factor(Compartments)))
	legend("topright",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi,DIM="12",indS,1)

	plot(PCAi$x[,1:2],cex=2,xlab=percenti[1],ylab=percenti[2],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_AFNP2)
	text(PCAi$x[,1],PCAi$x[,2],
		colnames(data.163)[indS],
		cex=1/2)
	legend("topleft",fill=hue_pal()(4),
		legend=levels(factor(Compartments)))
	legend("topright",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi,DIM="12",indS,2)

	plot(PCAi$x[,1:2],cex=2,xlab=percenti[1],ylab=percenti[2],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_levels2)
	legend("topleft",fill=seq(3),
		legend=levels(factor(Levels)))
	legend("topright",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi,DIM="12",indS,3)

#################################
#################################

	plot(PCAi$x[,1],PCAi$x[,3],cex=2,xlab=percenti[1],ylab=percenti[3],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_AFNP2)
	legend("topleft",fill=hue_pal()(4),
		legend=levels(factor(Compartments)))
	legend("bottomleft",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi,DIM="13",indS,1)


	plot(PCAi$x[,1],PCAi$x[,3],cex=2,xlab=percenti[1],ylab=percenti[3],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_AFNP2)
	text(PCAi$x[,1],PCAi$x[,3],
		colnames(data.163)[indS],
		cex=1/2)
	legend("topleft",fill=hue_pal()(4),
		legend=levels(factor(Compartments)))
	legend("bottomleft",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi, DIM="13",indS,2)


	plot(PCAi$x[,1],PCAi$x[,3],cex=2,xlab=percenti[1],ylab=percenti[3],
		main=MAINSTR,
		pch=PCH_ages2,
		col=COL_levels2)
	legend("topleft",fill=seq(3),
		legend=levels(factor(Levels)))
	legend("bottomleft",pch=c(16,17),pt.cex=2,
		legend=levels(factor(Ages)))
	SVMPLOT(PCAi,DIM="13",indS,3)

###
	MAINSTR2<-paste0("PCA on genes  ",IMPUTE,"\n(Based on ",
		length(geneSymb3),
		" genes with at least ",NUMsample," valid values)")

	geneFam<-sapply(geneSymb3,function(x){
		which(sapply(LNABA,function(y){
				x%in%y
		}))[1]
	})
	geneFam<-names(LNABA)[geneFam]
	geneFam[is.na(geneFam)]<-"others"
	print(geneFam)
	COL_geneFam<-hue_pal()(8)[as.integer(factor(geneFam))]
	plot(PCAj$x[,1:2],cex=2,xlab=percentj[1],ylab=percentj[2],
		main=MAINSTR2,
		pch=16,
		col=COL_geneFam)
	legend("topleft",fill=hue_pal()(8),
		legend=levels(factor(geneFam)),cex=1/2)
	text(PCAj$x[,1],PCAj$x[,2],
		geneSymb3,
		cex=1/2)
}

