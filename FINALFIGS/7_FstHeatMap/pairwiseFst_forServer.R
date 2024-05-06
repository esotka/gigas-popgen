library(hierfstat)
library(ape)
library(lattice)
library(RColorBrewer)
library(scales)

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

### Heatmaps with all pop'ns
rm(list=ls())
load("FINALFIGS/7_FstHeatMap/pairwiseWCfst.Rda") # name = wc
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
meta <- meta[!meta$pop%in%c("YOJ","MAI"),]
wc2 <- wc[rownames(wc)%in%meta$pop,rownames(wc)%in%meta$pop]

# sort by region
meta$Region2[meta$Region2=="nonSource"] = "Korea"
meta$Region2[meta$Region2=="Akkeshi"] = "Hokkaido"
meta$Region2[meta$Region2=="PNW"] = "noAmerica (WA/CA)"
meta$Region2[meta$Region2=="soCalifornia"] = "noAmerica (Cali)"
meta$Region2 = factor(meta$Region2)
meta$Region2 = factor(meta$Region2,levels=levels(meta$Region2)[c(4,10,12,5,3,11,8,7,6,9,1,2)])
siteorder = meta$pop[order(meta$Region2)]

wc3 = wc2[order(match(rownames(wc2),siteorder)),order(match(colnames(wc2),siteorder))]
colnames(wc3) <- paste(colnames(wc3),sort(meta$Region2))
#col.l <- colorRampPalette(brewer.pal(11, "RdBu")[11:1])
col.l = alpha(bluered(100),.5)
pdf('FINALFIGS/7_FstHeatMap/FstFigsTable_FINAL.pdf',width=15,height=13)
f <- levelplot(wc3,col.regions=col.l,,xlab="",ylab="",cuts=50)
print(f)
dev.off()
