# loci with PCA
#library(gplots)
library(hierfstat)
library(colorRamps)
library(scales)

rm(list=ls())
geno <- read.delim("FINALFIGS/SNPs_noZeros.txt",header=T)
rownames(geno) <- geno$id
geno <- geno[,-1]
# dim(geno)
#[1]  726 7045
pca1 <- prcomp(geno)
pdf("FINALFIGS/3_PCA/PC_allsamps.pdf")
pop <- substr(rownames(pca1$x),1,3)
out <- summary(pca1)
perc.pc1 <- round(out$importance[2,1],3)*100
perc.pc2 <- round(out$importance[2,2],3)*100
xbar <- aggregate(pca1$x[,1:2],by=list(pop),mean)
popRegPts <- xbar[match(pop,xbar$Group.1),]

col.to.use = c(blue2red(5),"darkgrey",rep("black",6))
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
popRegPts$reg <- factor(meta$Region2[match(popRegPts$Group.1,meta$pop)])
popRegPts$reg = factor(popRegPts$reg,levels(popRegPts$reg)[c(1,5,13,10,4,8,2,3,6,7,9,11,12)])
xbar$reg <- factor(meta$Region2[match(xbar$Group.1,meta$pop)])
xbar$reg = factor(xbar$reg,levels(xbar$reg)[c(1,5,13,10,4,8,2,3,6,7,9,11,12)])

plot(pca1$x[,1],pca1$x[,2],cex=0,xaxt="none",yaxt="none",ylab="",xlab="")
mtext(side=1,line=1,paste("PC1 (",perc.pc1,"%)",sep="")); mtext(side=2,line=1,paste("PC2 (",perc.pc2,"%)",sep=""))
segments(x0 = pca1$x[,1],y0=pca1$x[,2],x1=popRegPts$PC1,y1=popRegPts$PC2,col=alpha(col.to.use[popRegPts$reg],0.3))#col="lightgrey")#,col=alpha(col.sub,.5))
points(x=xbar$PC1,y=xbar$PC2,pch=20,cex=4,col=col.to.use[xbar$reg])
text(x=xbar$PC1,y=xbar$PC2,xbar$Group.1,cex=.5,col="white")
text(x=c(-6.2,-8.8,-6.4),y=c(-5.8,-0.4,2.3),c("Northern Europe","Argentina","Chile"),pos=4,col="red")
dev.off()

