# pretty HS

rm(list=ls())
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
stats <- read.csv("FINALFIGS/2_basicStats/BasicStats_all.csv")
pca_results <- read.csv("FINALFIGS/3_PCA/PC.out.csv")

pdf(width=9,height=4,"FINALFIGS/4_PC1&HS_vs_SST/PC1&Hs_prettyPlots.pdf")
#quartz(width=9,height=6)
par(mfrow=c(1,2),mar=c(4,5,1,1))

#### pc1 ~ sstmean
pca_results$sstmean <- meta$BO_sstmean[match(pca_results$Group.1,meta$pop)]
pca_results$natnon <- meta$NatNon[match(pca_results$Group.1,meta$pop)]
plot(PC1~sstmean,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],xlim=range(pca_results$sstmean),ylim=range(pca_results$PC1),col=pca_results$pc1.cols,pch=c(21,20)[factor(pca_results$natnon)],xlab="",ylab="")
abline(lm(PC1~sstmean,data=pca_results[pca_results$natnon=="Native",]))
abline(lm(PC1~sstmean,data=pca_results[pca_results$natnon=="Introduced",]),lty="dashed")

mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=pca_results$sstmean,y=pca_results$PC1,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)


### PC1 ~ Hs
pca_results$Hs <- stats$Hs[match(pca_results$Group.1,stats$X)]
pca_results2 <- pca_results[!pca_results$Group.1%in%c("MAI","YOJ"),] #n<5 for Hs

plot(PC1~Hs,data=pca_results2,cex=c(3,4)[factor(pca_results2$natnon)],xlim=range(pca_results2$Hs),ylim=range(pca_results2$PC1),col=pca_results2$pc1.cols,pch=c(21,20)[factor(pca_results2$natnon)],xlab="",ylab="")
abline(lm(PC1~Hs,data=pca_results2[pca_results2$natnon=="Native",]))
abline(lm(PC1~Hs,data=pca_results2[pca_results2$natnon=="Introduced",]),lty="dashed")

mtext("Expected Heterozygosity",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=pca_results$Hs,y=pca_results$PC1,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)

print(cor.test(~PC1+Hs,data=pca_results2[pca_results2$natnon=="Native",]))
print(cor.test(~PC1+Hs,data=pca_results2[pca_results2$natnon=="Introduced",]))

#### Hs ~ sstmean

plot(Hs~sstmean,data=pca_results2,cex=c(3,4)[factor(pca_results2$natnon)],xlim=range(pca_results2$sstmean),ylim=range(pca_results2$Hs,na.rm=T),col=pca_results2$pc1.cols,pch=c(21,20)[factor(pca_results2$natnon)],xlab="",ylab="")
abline(lm(Hs~sstmean,data=pca_results2[pca_results2$natnon=="Native",]))
abline(lm(Hs~sstmean,data=pca_results2[pca_results2$natnon=="Introduced",]),lty="dashed")

mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Expected Heterozygosity",side=2,line=2.5,cex=1.5)
text(x=pca_results$sstmean,y=pca_results$Hs,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)

### MAKE A BLANK

plot(Hs~sstmean,data=pca_results2,cex=c(3,4)[factor(pca_results2$natnon)],xlim=range(pca_results2$sstmean),ylim=range(pca_results2$Hs,na.rm=T),col=pca_results2$pc1.cols,pch=c(21,20)[factor(pca_results2$natnon)],xlab="",ylab="",type="n")

##########################
### Native only plots ####
##########################

nat = pca_results[pca_results$natnon=="Native",]

plot(PC1~sstmean,data=nat,cex=4,xlim=range(nat$sstmean),ylim=range(nat$PC1),col=nat$pc1.cols,pch=20,xlab="",ylab="")
abline(lm(PC1~sstmean,data=nat))

mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=nat$sstmean,y=nat$PC1,nat$Group.1,cex=.5,col=nat$popTextCol)


### PC1 ~ Hs 
nat2 <- nat[!nat$Group.1%in%c("MAI","YOJ"),] #n<5 for Hs

plot(PC1~Hs,data=nat2,cex=4,xlim=range(nat2$Hs),ylim=range(nat2$PC1),col=nat2$pc1.cols,pch=20,xlab="",ylab="")
abline(lm(PC1~Hs,data=nat2))

mtext("Expected Heterozygosity",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=nat$Hs,y=nat$PC1,nat$Group.1,cex=.5,col=nat$popTextCol)

##########################
### PC1 ~ Hs # no putative aquacultured pops (i.e., Argentina, Chile, northern Europe)
##########################
pca_results3 <- pca_results2[!pca_results2$GeneticRegions%in%c("Argentina","Chile","noEurope"),]

plot(PC1~Hs,data=pca_results3,cex=c(3,4)[factor(pca_results3$natnon)],xlim=range(pca_results3$Hs),ylim=range(pca_results3$PC1),col=pca_results3$pc1.cols,pch=c(21,20)[factor(pca_results3$natnon)],xlab="",ylab="")
abline(lm(PC1~Hs,data=pca_results3[pca_results3$natnon=="Native",]))
abline(lm(PC1~Hs,data=pca_results3[pca_results3$natnon=="Introduced",]),lty="dashed")

mtext("Expected Heterozygosity",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=pca_results3$Hs,y=pca_results3$PC1,pca_results3$Group.1,cex=.5,col=pca_results3$popTextCol)

print(cor.test(~PC1+Hs,data=pca_results3[pca_results3$natnon=="Native",]))
print(cor.test(~PC1+Hs,data=pca_results3[pca_results3$natnon=="Introduced",]))

##########################
### Hs ~ sstmean # no putative aquacultured pops (i.e., Argentina, Chile, northern Europe)
##########################
plot(Hs~sstmean,data=pca_results3,cex=c(3,4)[factor(pca_results3$natnon)],ylim=range(pca_results3$Hs),xlim=range(pca_results3$sstmean),col=pca_results3$pc1.cols,pch=c(21,20)[factor(pca_results3$natnon)],xlab="",ylab="")
abline(lm(Hs~sstmean,data=pca_results3[pca_results3$natnon=="Native",]))
abline(lm(Hs~sstmean,data=pca_results3[pca_results3$natnon=="Introduced",]),lty="dashed")

mtext("Expected Heterozygosity",side=2,line=2.5,cex=1.5)
mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
text(x=pca_results3$sstmean,y=pca_results3$Hs,pca_results3$Group.1,cex=.5,col=pca_results3$popTextCol)

print(cor.test(~Hs+sstmean,data=pca_results3[pca_results3$natnon=="Native",]))
print(cor.test(~Hs+sstmean,data=pca_results3[pca_results3$natnon=="Introduced",]))


dev.off()

