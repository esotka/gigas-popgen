library(strataG)
library(adegenet)
dat =read.table("dipsmll.351ind.geno_noMissingLoci.txt",header=F)
meta=read.csv("meta-SNP.pops_edited.V2.csv")

rownames(dat)=dat[,1]
dat=dat[,c(-1)]
pops=sapply(strsplit(rownames(dat),"_"),function(x)x[1])
dat=dat[pops%in%meta$pop,]

dat=dat[,sample(1:ncol(dat),1000,replace=F)] #select 1000 loci at random to populate gtype
calls = do.call(cbind,lapply(1:ncol(dat),function(i){ifelse(dat[,i]==0,"TT",ifelse(dat[,i]==1,"AT","AA"))}))
rownames(calls)=rownames(dat)
gin = genind2gtypes(df2genind(calls,ploidy=2,ncode=1))
strats = sapply(strsplit(rownames(calls),"_"),function(x){x[1]})
names(strats)=rownames(calls)
setStrata(gin) = strats

saveRDS(file="gverm_snps_229ind_1000loci.RDS",gin)
