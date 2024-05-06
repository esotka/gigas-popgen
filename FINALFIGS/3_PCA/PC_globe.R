### points of PC1 output of SNPs
rm(list=ls())
library(maps);library(fields);library(reshape);library(mapdata);library(colorRamps);library(plotrix);library(scales)
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
nat.reg <- factor(meta$GeneticRegions[match(nat.pop,meta$pop)])
nat.reg <- factor(nat.reg,levels(nat.reg)[c(1,4,6,5,2,3)])
non <- geno[natnon=="Introduced",]
non.pop <- substr(rownames(non),1,3)

pca1 <- prcomp(nat)
#quartz(width=6,height=6)
png("FINALFIGS/3_PCA/PC_globe.png",width=6,height=6,units="in",res=500)
par(mar=c(2,4,2,4))
#par(mfrow=c(1,2),mar=c(0,0,0,0))
plot(x=pca1$x[,2],y=pca1$x[,1],cex=0,xaxt="none",yaxt="none",ylab="",xlab="")
xbar <- aggregate(pca1$x[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$Region2[match(xbar$Group.1,meta$pop)]
#cols.to.use <- c(blue2red(5),"black")
#names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima","nonSource")
#xbar$pc1.cols = cols.to.use[match(xbar$GeneticRegions,names(cols.to.use))]
xbar$pc1.Rounded <- round(xbar$PC1+abs(min(xbar$PC1)),0)+1
#xbar$pc1.cols <- blue2red(n = max(xbar$pc1.Rounded))[max(xbar$pc1.Rounded):1][xbar$pc1.Rounded]
library(RColorBrewer)
#xbar$pc1.cols <- blue2red(n = max(xbar$pc1.Rounded))[max(xbar$pc1.Rounded):1][xbar$pc1.Rounded]
xbar$pc1.Rounded[xbar$pc1.Rounded>9] = 9
xbar$pc1.cols <- brewer.pal(n = 9,"Oranges")[1:max(xbar$pc1.Rounded)][xbar$pc1.Rounded]

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
xbar$popTextCol <- c("black","white")[factor(xbar$GeneticRegions=="nonSource")]
xbar$popTextCol[xbar$GeneticRegions=="Akkeshi"] <- "white"
text(x=-xbar$PC2,y=xbar$PC1,xbar$Group.1,cex=.5,col=xbar$popTextCol)
#popPts <- xbar.IntroPredict[match(non.pop,xbar.IntroPredict$Group.1),1:3]
#segments(x0 = -introPredict4[,2],y0=introPredict4[,1],x1=-popPts$PC2,y1=popPts$PC1,col=alpha(xbar.IntroPredict$pc1.cols,.1))
#points(x=-xbar.IntroPredict$PC2,y=xbar.IntroPredict$PC1,col=xbar.IntroPredict$pc1.cols,pch=21,lwd=2)

out <- summary(pca1)
perc.pc1 <- round(out$importance[2,1],3)*100
perc.pc2 <- round(out$importance[2,2],3)*100
mtext(side=2,line=-1,paste("PC1 (",perc.pc1,"%)",sep="")); mtext(side=1,line=-1,paste("PC2 (",perc.pc2,"%)",sep=""))
dev.off()

### maps
all <- rbind(xbar[,c("Group.1","PC1","PC2","pc1.cols")],xbar.IntroPredict[,c("Group.1","PC1","PC2","pc1.cols")])
all$Longitude <- meta$Longitude[match(all$Group.1,meta$pop)]
all$Latitude <- meta$Latitude[match(all$Group.1,meta$pop)]
all$lon.pie <- meta$lon.pie[match(all$Group.1,meta$pop)]
all$lat.pie <- meta$lat.pie[match(all$Group.1,meta$pop)]
all$GeneticRegions <- meta$Region2[match(all$Group.1,meta$pop)]
all$popTextCol <- c("black","white")[factor(all$GeneticRegions=="Korea")]
all$popTextCol[all$Group.1=="SAR"] <- "white"
# Japan
png('FINALFIGS/3_PCA/PC_globe_japan.png',width=5,height=4,units="in",res=500)
#quartz(width=5,height=4)
par(mar=c(0,0,0,0))
map("worldHires","China",bg="lightsteelblue1",fill=T,col="white",xlim=c(123,150),ylim=c(28,47),xaxt="n",yaxt="n",ylab="",xlab="")
map("worldHires","South Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","North Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","Japan",add=TRUE,col="white",fill=TRUE)
map("worldHires","USSR",col="white",fill=TRUE,add=TRUE)
box()
segments(all$Longitude,all$Latitude,all$lon.pie,all$lat.pie,col="red",lwd=2)
points(x=all$lon.pie,y=all$lat.pie,col=all$pc1.cols,pch=20,cex=4)
points(x=all$lon.pie,y=all$lat.pie,col="black",pch=21,cex=3)
text(x=all$lon.pie,y=all$lat.pie,all$Group.1,col=all$popTextCol,cex=.5)
dev.off()

### Europe range
png('FINALFIGS/3_PCA/PC_globe_europe.png',width=5,height=4,units="in",res=500)
#jpeg('Plot-points.in.color.non-native.jpg',width=10,height=8,units="in",res=500)
par(mar=c(0,0,0,0))
map('worldHires',xlim=c(-17,24),ylim=c(34,67),bg="lightsteelblue1",fill=T,col="white")
box()
segments(all$Longitude,all$Latitude,all$lon.pie,all$lat.pie,col="red",lwd=2)
points(x=all$lon.pie,y=all$lat.pie,col=all$pc1.cols,pch=20,cex=4)
points(x=all$lon.pie,y=all$lat.pie,col="black",pch=21,cex=3)
text(x=all$lon.pie,y=all$lat.pie,all$Group.1,col="black",cex=.5)
dev.off()

### globe
png('FINALFIGS/3_PCA/PC_globe_all.png',width=10,height=8,units="in",res=500)
map('worldHires')
tmp <- all[all$GeneticRegions%in%c("Argentina","Chile","New Zealand","PNW","soCalifornia"),]
lat <- tapply(tmp$Latitude,tmp$GeneticRegions,mean)
lon <- tapply(tmp$Longitude,tmp$GeneticRegions,mean)
cols.to.use <- tmp$pc1.cols[match(names(lat),tmp$GeneticRegions)]
points(x=lon,y=lat,col=cols.to.use,pch=20,cex=4)
points(x=lon,y=lat,col="black",pch=21,cex=3)
dev.off()

write.csv(all,"FINALFIGS/3_PCA/PC.out.csv",quote=F,row.names = F)
