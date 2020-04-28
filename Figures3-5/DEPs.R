agegrp<-"young"
OUTDIR<-"DEGs-for-young-DEC5"
outPDF<-"DEG.young.only.DEC4.pdf"
source("step5.DEGs.within.young.or.old.only.R")

##########################
agegrp<-"old"
OUTDIR<-"DEGs-for-old-DEC5"
outPDF<-"DEG.old.only.DEC4.pdf"
source("step5.DEGs.within.young.or.old.only.R")

##########################

source("step7.DEGs.young.vs.old.R")

res1<-fisher.test(matrix(c(6,3,19950-55-6,52),ncol=2))

res1<-fisher.test(matrix(c(44,55,19950-100-44,45),ncol=2))
