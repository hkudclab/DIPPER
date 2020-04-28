source("step5.DEGs.within.young.only.R")
source("step6.DEGs.within.old.only.R")
source("step7.DEGs.young.vs.old.R")

res1<-fisher.test(matrix(c(6,3,19950-55-6,52),ncol=2))

res1<-fisher.test(matrix(c(44,55,19950-100-44,45),ncol=2))
