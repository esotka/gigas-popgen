### points of PC1 output of SNPs
rm(list=ls())
library(maps);library(fields);library(reshape);library(mapdata);library(colorRamps);library(plotrix);library(scales);library(ggplot2);library(ggridges)
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
### PCA ####
geno <- read.delim("FINALFIGS/SNPs_noZeros.txt",header=T)
rownames(geno) <- geno$id
geno <- geno[,-1]
pop <- substr(rownames(geno),1,3)
natnon <- meta$NatNon[match(pop,meta$pop)]
#natnon <- ifelse(reg%in%c("Japan","Korea"),"Native","Introduced")
nat <- geno[natnon=="Native",]
nat.pop <- factor(substr(rownames(nat),1,3))
nat.reg <- factor(meta$Region2[match(nat.pop,meta$pop)])
nat.reg <- factor(nat.reg,levels(nat.reg)[c(1,3,6,5,2,4)])
non <- geno[natnon=="Introduced",]
non.pop <- substr(rownames(non),1,3)
non.reg <- meta$Region2[match(non.pop,meta$pop)]
pca1 <- prcomp(nat)
### predict ###
introPredict4 <- predict(pca1,newdata=non)

natnon <- data.frame(pop=c(as.character(nat.pop),non.pop),reg=c(as.character(nat.reg),non.reg),pc1=c(pca1$x[,1],introPredict4[,1]))
## remove Kag and non-source
natnon <- natnon[!(natnon$pop=="YOJ" | natnon$reg=="nonSource"),]
natnon$reg <- factor(natnon$reg)
natnon$reg <- factor(natnon$reg,levels=levels(natnon$reg)[c(1,4,11,8,5,3,2,7,9,6,10)])
## PRINT THIS TO A FILE
ggplot(natnon,aes(x = pc1, y = reg)) + geom_density_ridges_gradient(scale=.8) + coord_flip() 

#quartz(width=6,height=6)
png("FINALFIGS/3_PCA/PC_globe_separateNatNon.png",width=6,height=6,units="in",res=500)
par(mar=c(2,4,2,4))
#par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(x=pca1$x[,2],y=pca1$x[,1],cex=0,xaxt="none",yaxt="none",ylab="",xlab="")
xbar <- aggregate(pca1$x[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$GeneticRegions[match(xbar$Group.1,meta$pop)]
xbar$pc1.Rounded <- round(xbar$PC1+abs(min(xbar$PC1)),0)+1
xbar$pc1.cols <- blue2red(n = max(xbar$pc1.Rounded))[max(xbar$pc1.Rounded):1][xbar$pc1.Rounded]
### predict ###
introPredict4 <- predict(pca1,newdata=non)
xbar.IntroPredict <- aggregate(introPredict4[,1:2],by=list(non.pop),mean)
xbar.IntroPredict$pc1.Rounded <- round(xbar.IntroPredict$PC1+abs(min(xbar$PC1)),0)+1
xbar.IntroPredict$pc1.cols <- xbar$pc1.cols[match(xbar.IntroPredict$pc1.Rounded,xbar$pc1.Rounded)]

xbar$pc1.cols[xbar$Group.1=="MAI"] <- "grey"
xbar$pc1.cols[xbar$GeneticRegions=="Korea"] <- "black"

popPts <- xbar[match(nat.pop,xbar$Group.1),]
seg_col <- xbar$pc1.cols[match(nat.pop,as.character(xbar$Group.1))]
segments(x0 = -pca1$x[,2],y0=pca1$x[,1],x1=-popPts$PC2,y1=popPts$PC1,col=alpha(seg_col,0.2))#,col=alpha(col.sub,.5))
points(x=-xbar$PC2,y=xbar$PC1,pch=20,cex=4,col=xbar$pc1.cols)
xbar$popTextCol <- c("black","white")[factor(xbar$GeneticRegions=="Korea")]
xbar$popTextCol[xbar$Group.1=="SAR"] <- "white"
text(x=-xbar$PC2,y=xbar$PC1,xbar$Group.1,cex=.5,col=xbar$popTextCol)
popPts <- xbar.IntroPredict[match(non.pop,xbar.IntroPredict$Group.1),1:3]
segments(x0 = -introPredict4[,2],y0=introPredict4[,1],x1=-popPts$PC2,y1=popPts$PC1,col=alpha(xbar.IntroPredict$pc1.cols,.1))
points(x=-xbar.IntroPredict$PC2,y=xbar.IntroPredict$PC1,col=xbar.IntroPredict$pc1.cols,pch=21,lwd=2)

out <- summary(pca1)
perc.pc1 <- round(out$importance[2,1],3)*100
perc.pc2 <- round(out$importance[2,2],3)*100
mtext(side=2,line=-1,paste("PC1 (",perc.pc1,"%)",sep="")); mtext(side=1,line=-1,paste("PC2 (",perc.pc2,"%)",sep=""))
dev.off()

