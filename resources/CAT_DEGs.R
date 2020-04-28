load("G:/kathy-cheah-ddd/resources/GO_INFLAMMATORY_RESPONSE.RData")
load("G:/kathy-cheah-ddd/resources/ligands.GO.and.HGNC.RData")
load("G:/kathy-cheah-ddd/resources/NABA.break.down.RData")
load("G:/kathy-cheah-ddd/resources/TF.CD.receptors.human.mouse.RData")

library(grid)
library(lattice)
library(gridExtra)
library(stringr)


CAT_DEGs<-function(MATin, BSIZE=5, WIDTH=80,Xpos=0.75){
		COLLAGENS<-as.character(MATin$FeatureName)[grep("^COL[0-9]+",
			as.character(MATin$FeatureName),
			ignore.case = T)]
		genei<- as.character(MATin$FeatureName)
		geneCat<-function(g1,g2){
			g3<-g1[toupper(g1)%in%toupper(as.character(g2))]
			paste0(g3,collapse=", ")
		}

##################
		gD<-unique(genei[indSig.down]);
		gU<-unique(genei[indSig.up]);
		LD<-length(gD)
		LU<-length(gU)
		WD<-WIDTH*2*(LD/(LD+LU))
		WU<-WIDTH*2*(LU/(LD+LU))

##################
		TFdown<-geneCat(gD,TF.human)
		CDdown<-geneCat(gD,surface.human)

		INFLdown<-geneCat(gD,GO_INFLAMMATORY_RESPONSE)
		LIGANDSdown<-geneCat(gD,c(ligands.GO,ligands.HGNC1,ligands.HGNC2))
		COLLAGENSdown<-geneCat(gD,COLLAGENS)
		PROTEOGLYCANSdown<-geneCat(gD,LNABA$NABA_PROTEOGLYCANS)
		GLYCOPROTEINSdown<-geneCat(gD,LNABA$NABA_ECM_GLYCOPROTEINS)
		ECM_REGULATORSdown<-geneCat(gD,LNABA$NABA_ECM_REGULATORS)
		ECM_AFFILIATEDdown<-geneCat(gD,LNABA$NABA_ECM_AFFILIATED)
		SECRETED_FACTORSdown<-geneCat(gD,LNABA$NABA_SECRETED_FACTORS)

		str1<-c(TFdown,CDdown,INFLdown,LIGANDSdown,COLLAGENSdown,PROTEOGLYCANSdown,
			GLYCOPROTEINSdown,ECM_REGULATORSdown,ECM_AFFILIATEDdown,SECRETED_FACTORSdown)

		othergene<-setdiff(gD,unlist(strsplit(str1,", ")))
		otherdown<-paste0(othergene,collapse=", ")
		##################

		TFup<-geneCat(gU,TF.human)
		CDup<-geneCat(gU,surface.human)

		INFLup<-geneCat(gU,GO_INFLAMMATORY_RESPONSE)
		LIGANDSup<-geneCat(gU,c(ligands.GO,ligands.HGNC1,ligands.HGNC2))
		COLLAGENSup<-geneCat(gU,COLLAGENS)
		PROTEOGLYCANSup<-geneCat(gU,LNABA$NABA_PROTEOGLYCANS)
		GLYCOPROTEINSup<-geneCat(gU,LNABA$NABA_ECM_GLYCOPROTEINS)
		ECM_REGULATORSup<-geneCat(gU,LNABA$NABA_ECM_REGULATORS)
		ECM_AFFILIATEDup<-geneCat(gU,LNABA$NABA_ECM_AFFILIATED)
		SECRETED_FACTORSup<-geneCat(gU,LNABA$NABA_SECRETED_FACTORS)

		str1<-c(TFup,CDup,INFLup,LIGANDSup,COLLAGENSup,PROTEOGLYCANSup,
			GLYCOPROTEINSup,ECM_REGULATORSup,ECM_AFFILIATEDup,SECRETED_FACTORSup)

		othergene<-setdiff(gU,unlist(strsplit(str1,", ")))
		otherup<-paste0(othergene,collapse=", ")
		##################
		WD<-WIDTH/4+WIDTH*1.5*(LD/(LD+LU))
		WU<-WIDTH/4+WIDTH*1.5*(LU/(LD+LU))
		##################
		mytheme <- gridExtra::ttheme_default(base_size=BSIZE,
			padding = unit(c(1, 3), "mm"))

		mat1<-rbind(c(str_wrap(TFdown,WD),str_wrap(TFup,WU)),
			c(str_wrap(CDdown,WD),str_wrap(CDup,WU)),
			c(str_wrap(COLLAGENSdown,WD),str_wrap(COLLAGENSup,WU)),
			c(str_wrap(PROTEOGLYCANSdown,WD),str_wrap(PROTEOGLYCANSup,WU)),
			c(str_wrap(GLYCOPROTEINSdown,WD),str_wrap(GLYCOPROTEINSup,WU)),
			c(str_wrap(ECM_REGULATORSdown,WD),str_wrap(ECM_REGULATORSup,WU)),
			c(str_wrap(ECM_AFFILIATEDdown,WD),str_wrap(ECM_AFFILIATEDup,WU)),
			c(str_wrap(SECRETED_FACTORSdown,WD),str_wrap(SECRETED_FACTORSup,WU)),
			c(str_wrap(LIGANDSdown,WD),str_wrap(LIGANDSup,WU)),
			c(str_wrap(INFLdown,WD),str_wrap(INFLup,WU)),
			c(str_wrap(otherdown,WD),str_wrap(otherup,WU)))
		rownames(mat1)<-c("TF","CD","Collagens",
			"Proteoglycans",
			"Glycoproteins",
			"ECM Regulators",
			"ECM Affiliated",
			"Secreted factors",
			"Ligands",
			"Inflammatory",
			"others")
		colnames(mat1)<-c("up in blue","up in red")
		pushViewport(viewport(x=Xpos,y=0.5))
		grid.table(mat1,theme=mytheme )
		popViewport()
}



CAT_DEG2<-function(gUP=c(""),gDown=c(""), BSIZE=5, WIDTH=80,Xpos=0.75){

		gD<-unique(gDown);
		gU<-unique(gUP);
		LD<-length(gD)
		LU<-length(gU)
		WD<-WIDTH*2*(LD/(LD+LU))
		WU<-WIDTH*2*(LU/(LD+LU))

		gIn<-c(gUP,gDown)
		COLLAGENS<-gIn[grep("^COL[0-9]+",gIn,ignore.case = T)]

		geneCat<-function(g1,g2){
			g3<-g1[toupper(g1)%in%toupper(as.character(g2))]
			g3<-setdiff(unique(g3),"")
			paste0(g3,collapse=", ")
		}

##################
		gi<-gD;
		TFdown<-geneCat(gi,TF.human)
		CDdown<-geneCat(gi,surface.human)

		INFLdown<-geneCat(gi,GO_INFLAMMATORY_RESPONSE)
		LIGANDSdown<-geneCat(gi,c(ligands.GO,ligands.HGNC1,ligands.HGNC2))
		COLLAGENSdown<-geneCat(gi,COLLAGENS)
		PROTEOGLYCANSdown<-geneCat(gi,LNABA$NABA_PROTEOGLYCANS)
		GLYCOPROTEINSdown<-geneCat(gi,LNABA$NABA_ECM_GLYCOPROTEINS)
		ECM_REGULATORSdown<-geneCat(gi,LNABA$NABA_ECM_REGULATORS)
		ECM_AFFILIATEDdown<-geneCat(gi,LNABA$NABA_ECM_AFFILIATED)
		SECRETED_FACTORSdown<-geneCat(gi,LNABA$NABA_SECRETED_FACTORS)

		str1<-c(TFdown,CDdown,INFLdown,LIGANDSdown,COLLAGENSdown,PROTEOGLYCANSdown,
			GLYCOPROTEINSdown,ECM_REGULATORSdown,ECM_AFFILIATEDdown,SECRETED_FACTORSdown)

		othergene<-setdiff(gi,unlist(strsplit(str1,", ")))
		otherdown<-paste0(othergene,collapse=", ")
		##################
		gi<-gU

		TFup<-geneCat(gi,TF.human)
		CDup<-geneCat(gi,surface.human)

		INFLup<-geneCat(gi,GO_INFLAMMATORY_RESPONSE)
		LIGANDSup<-geneCat(gi,c(ligands.GO,ligands.HGNC1,ligands.HGNC2))
		COLLAGENSup<-geneCat(gi,COLLAGENS)
		PROTEOGLYCANSup<-geneCat(gi,LNABA$NABA_PROTEOGLYCANS)
		GLYCOPROTEINSup<-geneCat(gi,LNABA$NABA_ECM_GLYCOPROTEINS)
		ECM_REGULATORSup<-geneCat(gi,LNABA$NABA_ECM_REGULATORS)
		ECM_AFFILIATEDup<-geneCat(gi,LNABA$NABA_ECM_AFFILIATED)
		SECRETED_FACTORSup<-geneCat(gi,LNABA$NABA_SECRETED_FACTORS)

		str1<-c(TFup,CDup,INFLup,LIGANDSup,COLLAGENSup,PROTEOGLYCANSup,
			GLYCOPROTEINSup,ECM_REGULATORSup,ECM_AFFILIATEDup,SECRETED_FACTORSup)

		othergene<-setdiff(gi,unlist(strsplit(str1,", ")))
		otherup<-paste0(othergene,collapse=", ")
		##################
		mytheme <- gridExtra::ttheme_default(base_size=BSIZE,
			padding = unit(c(1, 3), "mm"))

		mat1<-rbind(c(str_wrap(TFdown,WD),str_wrap(TFup,WU)),
			c(str_wrap(CDdown,WD),str_wrap(CDup,WU)),
			c(str_wrap(COLLAGENSdown,WD),str_wrap(COLLAGENSup,WU)),
			c(str_wrap(PROTEOGLYCANSdown,WD),str_wrap(PROTEOGLYCANSup,WU)),
			c(str_wrap(GLYCOPROTEINSdown,WD),str_wrap(GLYCOPROTEINSup,WU)),
			c(str_wrap(ECM_REGULATORSdown,WD),str_wrap(ECM_REGULATORSup,WU)),
			c(str_wrap(ECM_AFFILIATEDdown,WD),str_wrap(ECM_AFFILIATEDup,WU)),
			c(str_wrap(SECRETED_FACTORSdown,WD),str_wrap(SECRETED_FACTORSup,WU)),
			c(str_wrap(LIGANDSdown,WD),str_wrap(LIGANDSup,WU)),
			c(str_wrap(INFLdown,WD),str_wrap(INFLup,WU)),
			c(str_wrap(otherdown,WD),str_wrap(otherup,WU)))
		rownames(mat1)<-c("TF","CD","Collagens",
			"Proteoglycans",
			"Glycoproteins",
			"ECM Regulators",
			"ECM Affiliated",
			"Secreted factors",
			"Ligands",
			"Inflammatory",
			"others")
		colnames(mat1)<-c(paste0("up in blue (n=",length(gDown),")"),
			paste0("up in red (n=",length(gUP),")"))
		pushViewport(viewport(x=Xpos,y=0.5))
		grid.table(mat1,theme=mytheme )
		popViewport()
}