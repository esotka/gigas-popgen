# pretty HS

rm(list=ls())
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
stats <- read.csv("FINALFIGS/2_basicStats/BasicStats_all.csv")
#pca_results <- read.csv("FINALFIGS/3_PCA/PC.out.csv")

### PCA - native pops ####
geno <- read.delim("FINALFIGS/SNPs_noZeros.txt",header=T)
rownames(geno) <- geno$id
geno <- geno[,-1]
pop <- substr(rownames(geno),1,3)
natnon <- meta$NatNon[match(pop,meta$pop)]
nat <- geno[natnon=="Native",]
nat.pop <- factor(substr(rownames(nat),1,3))
nat.reg <- factor(meta$Region2[match(nat.pop,meta$pop)])
nat.reg <- factor(nat.reg,levels(nat.reg)[c(1,3,6,5,2,4)])
pca1 <- prcomp(nat)
pca1_mat = pca1$x[,1:2]
xbar <- aggregate(pca1_mat[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$Region2[match(xbar$Group.1,meta$pop)]
cols.to.use <- c(blue2red(5),"black")
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima","nonSource")
xbar$cols = cols.to.use[match(xbar$GeneticRegions,names(cols.to.use))]
#xbar$popTextCol = "black"
#xbar$popTextCol[xbar$Group.1%in%c("SAR","GOS","GWG","KOJ","MAI","AKK")] <- "white"
### Assign PC1 - Non-native ###
pop_non <- meta$pop[meta$NatNon=="Introduced"]
non = geno[pop%in%pop_non,]
non.pop <- substr(rownames(non),1,3)
non.reg <- meta$Region2[match(non.pop,meta$pop)]
introPredict <- predict(pca1,newdata=non)[,1:2]

xbar_non <- aggregate(introPredict[,1:2],by=list(non.pop),mean)
xbar_non$GeneticRegions <- meta$Region2[match(xbar_non$Group.1,meta$pop)]
xbar_non$cols="black"
pca_results = rbind(xbar,xbar_non)
pca_results <- pca_results[!pca_results$Group.1%in%c("MAI","YOJ"),] #n<5


pdf(width=9,height=4,"FINALFIGS/4_PC1&HS_vs_SST/PC1&Hs_prettyPlots.pdf")
#quartz(width=9,height=6)
par(mfrow=c(1,2),mar=c(4,5,1,1))

#### pc1 ~ sstmean
pca_results$sstmean <- meta$BO_sstmean[match(pca_results$Group.1,meta$pop)]
pca_results$natnon <- meta$NatNon[match(pca_results$Group.1,meta$pop)]
pca_results$reg = meta$Region2[match(pca_results$Group.1,meta$pop)]
pca_results$reg2 = as.factor(pca_results$reg)
pca_results$symbols = c(20,21,21,20,23,21,20,21,20,22,21,20)[pca_results$reg2]
#pca_results$symbols[is.na(pca_results$symbols)] = 20
plot(PC1~sstmean,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],
    xlim=c(6,22),ylim=range(pca_results$PC1),
    col=pca_results$cols,pch=pca_results$symbols,xlab="",ylab="")
points(x=pca_results$sstmean[pca_results$natnon=="Native"],y=pca_results$PC1[pca_results$natnon=="Native"],pch=21,cex=2.9)
abline(lm(PC1~sstmean,data=pca_results[pca_results$natnon=="Native",]))
abline(lm(PC1~sstmean,data=pca_results[pca_results$natnon=="Introduced",]),lty="dashed")

mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
#text(x=pca_results$sstmean,y=pca_results$PC1,pca_results$Group.1,cex=.5)#,col=pca_results$popTextCol)

print(cor.test(~PC1+sstmean,data=pca_results[pca_results$natnon=="Native",]))
print(cor.test(~PC1+sstmean,data=pca_results[pca_results$natnon=="Introduced",]))

### PC1 ~ Hs
pca_results$Hs <- stats$Hs[match(pca_results$Group.1,stats$X)]

plot(PC1~Hs,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],
    xlim=range(pca_results$Hs),ylim=range(pca_results$PC1),
    col=pca_results$cols,pch=pca_results$symbols,xlab="",ylab="")
points(x=pca_results$Hs[pca_results$natnon=="Native"],y=pca_results$PC1[pca_results$natnon=="Native"],pch=21,cex=2.9)

abline(lm(PC1~Hs,data=pca_results[pca_results$natnon=="Native",]))
abline(lm(PC1~Hs,data=pca_results[pca_results$natnon=="Introduced",]),lty="dashed")

mtext("Expected Heterozygosity",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
#text(x=pca_results$Hs,y=pca_results$PC1,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)

print(cor.test(~PC1+Hs,data=pca_results[pca_results$natnon=="Native",]))
print(cor.test(~PC1+Hs,data=pca_results[pca_results$natnon=="Introduced",]))

#### Hs ~ sstmean

plot(Hs~sstmean,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],xlim=range(pca_results$sstmean),
    ylim=range(pca_results$Hs,na.rm=T),col=pca_results$cols,pch=c(21,20)[factor(pca_results$natnon)],xlab="",ylab="")
abline(lm(Hs~sstmean,data=pca_results[pca_results$natnon=="Native",]))
abline(lm(Hs~sstmean,data=pca_results[pca_results$natnon=="Introduced",]),lty="dashed")

mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Expected Heterozygosity",side=2,line=2.5,cex=1.5)
#text(x=pca_results$sstmean,y=pca_results$Hs,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)

print(cor.test(~Hs+sstmean,data=pca_results[pca_results$natnon=="Native",]))
print(cor.test(~Hs+sstmean,data=pca_results[pca_results$natnon=="Introduced",]))

dev.off()

