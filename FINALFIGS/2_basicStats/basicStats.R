## from hierfstat
rm(list=ls())
library(hierfstat)
library(ggplot2)
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

#####
out.summary <- read.csv("FINALFIGS/2_basicStats/BasicStats_all.csv")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
out.summary$NatNon <- meta$NatNon[match(out.summary$X,meta$pop)]
out.summary$region <- meta$Region2[match(out.summary$X,meta$pop)]
anova(lm(Hs~NatNon,out.summary))
anova(lm(Ho~NatNon,out.summary))
anova(lm(Fis~NatNon,out.summary))

### three regions

out.summary$NatNonAqua <- out.summary$NatNon
out.summary$NatNonAqua[out.summary$region%in%c("Argentina","Chile","noEurope")] = "Introduced (Arg, Chi, NEur)"
out.summary$NatNonAqua[out.summary$NatNonAqua=="Introduced"] = "Introduced (other)"
#anova(lm(Hs~NatNonAqua,out.summary))
#anova(lm(Ho~NatNonAqua,out.summary))
#anova(lm(Fis~NatNonAqua,out.summary))
m <- lm(Hs~NatNonAqua,out.summary)
print(anova(m))
print(TukeyHSD(aov(m)))
library(pgirmess)
print(kruskal.test(Hs~NatNonAqua,out.summary))
library(pgirmess)
print(kruskalmc(Hs~NatNonAqua,out.summary))

pdf("FINALFIGS/2_basicStats/Hs.pdf",width=5,height=5)


f5 <- ggplot(out.summary,aes(y=Hs,x=NatNonAqua)) +
  geom_boxplot() +
  #geom_point(pch=20,cex=.3) + 
  geom_jitter(width = .05, alpha = .5) +
  xlab("") + ylab("Expected Heterozygosity") + ylim(c(0.135,0.160)) +
  guides(fill="none") +
  theme_classic() +
  annotate("text", x=c(1.5,2,2.5), y=c(0.151,0.156,0.151), label= c("ns","*","ns")) +
  annotate("segment", x=1.1, xend=1.9, y=0.15, yend=0.15) +
  annotate("segment", x=1.1, xend=2.9, y=0.155, yend=0.155) +
  annotate("segment", x=2.1, xend=2.9, y=0.15, yend=0.15)

plot(f5)
dev.off()


