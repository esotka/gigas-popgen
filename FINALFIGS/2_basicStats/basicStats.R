## from hierfstat
rm(list=ls())
library(hierfstat)
snp <- read.delim("FINALFIGS/SNPs_noZeros.txt",sep="\t")
rownames(snp) <- snp[,1]; snp <- snp[,-1]
pop <- substr(rownames(snp),1,3)

geno.h <- data.frame(grp=pop,snp)
out <- basic.stats(geno.h) # takes a minute

n <- c(table(pop))
Ho <- colMeans(out$Ho,na.rm=T) # observed heterozygosities
Hs <- colMeans(out$Hs,na.rm=T) # observed gene diversities ("sometimes misleadingly called expected heterozygosity")
Fis <- colMeans(out$Fis,na.rm=T) # observed Fis ==> these were all NAs

out.summary <- data.frame(n,Ho,Hs,Fis)
write.csv(out.summary,"FINALFIGS/2_basicStats/BasicStats_all.csv")

