### LD and Ne
library(strataG)
library(ggplot2)
rm(list=ls())
#snp <- read.delim("FINALFIGS/SNPs_noZeros.txt",sep="\t")
### convert into 2 alleles per locus
#nloci <- ncol(snp)-1
############################
### do this code once ######
#snp2 <- c()
#for (i in 1:nloci)
#{
#  tmp <- snp[,1+i]
#  key <- data.frame(snp=c(0,1,2),a=c(1,1,2),b=c(1,2,2))
#  a <- key$a[match(tmp,key$snp)]
#  b <- key$b[match(tmp,key$snp)]
#  snp2 <- cbind(snp2,a,b)
#}
#colnames(snp2) <- sort(c(paste("a",1:nloci,sep=""),paste("b",1:nloci,sep="")))
#pop <- substr(snp$id,1,3)
#dat <- data.frame(snp$id,pop,snp2)
#colnames(dat) <- c("Ind","Pop",colnames(snp2))
#gi <- df2gtypes(dat, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
#saveRDS(gi,file="FINALFIGS/cgiga_snp_gtype") ### gtype-formatted data
gin = readRDS(file="FINALFIGS/cgiga_snp_gtype")

pops <- unique(getStrata(gin))
pops <- pops[!pops=="YOJ"]
out <- c()
for (i in 1:length(pops))
{
gin2 <- gin[getStrata(gin)==pops[i]]
#LDgenepop(gin2)
out <- rbind(out,data.frame(
  ldNe(gin2), 
  hs=mean(heterozygosity(gin2)$exptd.het),
  ho = mean(heterozygosity(gin2,type = "observed")[,2])))
}

n = table(getStrata(gin))
n = n[!names(n)=="YOJ"]
out2 <- data.frame(n,out)

meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
out2$lat <- meta$Latitude[match(out2$Var1,meta$pop)]
out2$NatNon <- meta$NatNon[match(out2$Var1,meta$pop)]
out2$region <- meta$Region[match(out2$Var1,meta$pop)]
out2$region2 <- meta$Region2[match(out2$Var1,meta$pop)]
write.csv(out2,"FINALFIGS/2_basicStats/LD_and_Ne.csv",quote=F)
#print(ggscatmat(data=out2[,c("Ne","hs","ho","lat")],alpha=.1)) # PC1 and biotic  variables

####
out2 <- read.csv("FINALFIGS/2_basicStats/LD_and_Ne.csv")

tmp <- out2[!is.infinite(out2$Ne),] # remove Ne that are infinite.
 # GOS - Korea
 # OHK - Tokyo
 # PES - so Cal
 # TJE - so Cal
 # WLB - PNW
pdf("FINALFIGS/2_basicStats/LD_and_Ne.pdf",width=5,height=5)
f5 <- ggplot(tmp,aes(y=Ne,x=NatNon)) +#,ymax=param.uci,ymin=param.lci)) +
  geom_boxplot() +
  geom_point(pch=20,cex=.3) + 
  geom_jitter(width = .05, alpha = .5) +
  xlab("") +
  guides(fill=FALSE) +
  theme_classic() 
plot(f5)
dev.off()
anova(lm(Ne~NatNon,tmp))


