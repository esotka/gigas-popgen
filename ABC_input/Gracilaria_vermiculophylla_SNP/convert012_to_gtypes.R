library(strataG)
dat =read.table("dipsmll.351ind.geno_noMissingLoci.txt",header=F)

rownames(dat)=dat[,1]
dat=dat[,-1]
lst = lapply(1:ncol(dat),function(i){ifelse(dat[,i]==0,"TT",ifelse(dat[,i]==1,"AT","AA"))})
