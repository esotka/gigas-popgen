### PCOA on genetic distance
rm(list=ls())
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
#meta <- meta[!meta$Region2%in%c("nonSource","Kagoshima"),]
meta <- meta[!meta$pop%in%c("YOJ","MAI"),]

snp <- read.delim("FINALFIGS/SNPs_noZeros.txt",sep="\t")
rownames(snp) <- snp[,1]; snp <- snp[,-1]
pop <- substr(rownames(snp),1,3)
chr <- unlist(lapply(strsplit(colnames(snp),"_"),"[[",2))
snp2 <- snp[pop%in%meta$pop,]
pop2 <- substr(rownames(snp2),1,3)
geno.h <- data.frame(grp=pop2,snp2)
nei_dist = genet.dist(geno.h,method = "Da") # nei genetic distance
pcoa = dudi.pco(nei_dist,scannf=F,nf=2) # PCoA

plot(pcoa$tab[,1:2],type="n")
text(pcoa$tab[,1], pcoa$tab[,2], as.character(attr(nei_dist,"Labels")))

#out.dd <- as.dist(nei_dist)
hc <- hclust(nei_dist,method="ward.D2") # make new hclust
plot(hc,xlab="",sub="",cex=1.5)
print(rect.hclust(hc,k=5))

print(rect.hclust(hc,k=5))
out <- unlist(rect.hclust(hc,k=5)) # order levelplot() by clusters
siteorder <- data.frame(pop=names(out),num=1:length(out))
siteorder.mat <- siteorder$num[match(rownames(as.matrix(nei_dist)),siteorder$pop)]
wc <- as.matrix(nei_dist)[order(siteorder.mat),order(siteorder.mat)]

levelplot(wc,cex=.1,xlab="",ylab="")
