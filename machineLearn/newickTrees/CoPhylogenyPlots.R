### cophylogeny plots
library(phytools)
library(ape)
rm(list=ls())
# gigas SNPs
tr.root.giga <- read.tree("gigas_snp.tr")
## gverm usat
tr.root.gverm <- read.tree("gverm_usat.tr")

pdf("cgigasSNP_vs_others.pdf",width=7,height=10)
par(mar=c(1,1,1,1),mfrow=c(3,2))
obj<-cophylo(tr.root.giga,tr.root.gverm,rotate=F)
plot(obj)
mtext(side=3,at=c(-.3,.3),line=-3,c("gigas SNPs","gverm usat"),col="red",cex=.7)

#dev.off()

#battr mtDNA - from dendrogram of PhiST
tr.battr <- read.tree("battr_mtDNA.tr")
tr.root.battr <- root(tr.battr,outgroup="nonSource")
#plot(tr.root.battr)
#pdf("cgigasSNP_battrMTDNA.pdf",width=5,height=5)
#par(mar=c(1,1,1,1))
obj<-cophylo(tr.root.giga,tr.root.battr,rotate=F)
plot(obj)
mtext(side=3,at=c(-.3,.3),line=-3,c("gigas SNPs","battr mtDNA"),col="red",cex=.7)
#dev.off()

#upinn mtDNA - from dendrogram of PhiST
tr.upinn <- read.tree("upinn_mtDNA.tr")
tr.root.upinn <- root(tr.upinn,outgroup="nonSource")
#plot(tr.root.upinn)
#pdf("cgigasSNP_upinnMTDNA.pdf",width=5,height=5)
#par(mar=c(1,1,1,1))
obj<-cophylo(tr.root.giga,tr.root.upinn,rotate=F)
plot(obj)
mtext(side=3,at=c(-.3,.3),line=-3,c("gigas SNPs","upinn mtDNA"),col="red",cex=.7)
#dev.off()

#smuti usat - from dendrogram of PhiST
tr.smuti <- read.tree("smuti_usat.tr")
tr.root.smuti <- root(tr.smuti,outgroup="Russia")
#plot(tr.root.smuti)
#pdf("cgigasSNP_smutiusatt.pdf",width=5,height=5)
#par(mar=c(1,1,1,1))
obj<-cophylo(tr.root.giga,tr.root.smuti,rotate=F)
plot(obj)
mtext(side=3,at=c(-.3,.3),line=-3,c("gigas SNPs","smuti usat"),col="red",cex=.7)
dev.off()
