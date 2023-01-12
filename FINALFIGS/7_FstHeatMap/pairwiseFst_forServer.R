library(hierfstat)
library(ape)
library(lattice)
library(RColorBrewer)

rm(list=ls())
snp <- read.delim("FINALFIGS/SNPs_noZeros.txt",sep="\t")
rownames(snp) <- snp[,1]; snp <- snp[,-1]
pop <- substr(rownames(snp),1,3)
chr <- unlist(lapply(strsplit(colnames(snp),"_"),"[[",2))
snp2 <- snp[!pop%in%c("YOJ","MAI"),]
pop2 <- substr(rownames(snp2),1,3)
geno.h <- data.frame(grp=pop2,snp2)
#wc <- pairwise.WCfst(geno.h) ### takes a long time
#save(wc,"FINALFIGS/7_FstHeatMap/pairwiseWCfst.Rda")


### Heatmaps without putative aquacultured pop'ns
load("FINALFIGS/7_FstHeatMap/pairwiseWCfst.Rda")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
toExclude <- c("Argentina","Chile","Denmark","Ireland","Sweden","Norway")
meta <- meta[!meta$Region%in%toExclude,]
meta <- meta[!meta$pop=="MAI",]
wc2 <- wc[rownames(wc)%in%meta$pop,rownames(wc)%in%meta$pop]

### lattice::levelplot
out.dd <- as.dist(wc2)
hc <- hclust(out.dd,method = "average") # make new hclust
pdf('FINALFIGS/7_FstHeatMap/FstFigsTable.pdf',width=14,height=12)
plot(hc,xlab="",sub="",cex=1.5)
print(rect.hclust(hc,k=5))
out <- unlist(rect.hclust(hc,k=5)) # order levelplot() by clusters
siteorder <- data.frame(pop=names(out),num=1:length(out))
siteorder.wc2 <- siteorder$num[match(rownames(wc2),siteorder$pop)]
wc4 <- wc2[order(siteorder.wc2),order(siteorder.wc2)]
reg.hc <- meta$Region[match(rownames(wc4),meta$pop)]
colnames(wc4) <- paste(colnames(wc4),reg.hc)
col.l <- colorRampPalette(brewer.pal(11, "RdBu")[11:1])
f <- levelplot(wc4,cex=.1,col.regions=col.l,xlab="",ylab="")
print(f)
dev.off()

### Heatmaps with all pop'ns
rm(list=ls())
load("FINALFIGS/7_FstHeatMap/pairwiseWCfst.Rda")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
meta <- meta[!meta$pop=="MAI",]
wc2 <- wc[rownames(wc)%in%meta$pop,rownames(wc)%in%meta$pop]

### lattice::levelplot
out.dd <- as.dist(wc2)
hc <- hclust(out.dd,method = "average") # make new hclust
pdf('FINALFIGS/7_FstHeatMap/FstFigsTable_allpops.pdf',width=15,height=13)
plot(hc,xlab="",sub="",cex=1.5)
print(rect.hclust(hc,k=5))
out <- unlist(rect.hclust(hc,k=5)) # order levelplot() by clusters
siteorder <- data.frame(pop=names(out),num=1:length(out))
siteorder.wc2 <- siteorder$num[match(rownames(wc2),siteorder$pop)]
wc4 <- wc2[order(siteorder.wc2),order(siteorder.wc2)]
reg.hc <- meta$Region[match(rownames(wc4),meta$pop)]
colnames(wc4) <- paste(colnames(wc4),reg.hc)
col.l <- colorRampPalette(brewer.pal(11, "RdBu")[11:1])
f <- levelplot(wc4,cex=.1,col.regions=col.l,xlab="",ylab="")
print(f)
dev.off()




