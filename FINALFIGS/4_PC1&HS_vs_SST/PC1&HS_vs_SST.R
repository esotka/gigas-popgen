# pretty HS

rm(list=ls())
#library(mapdata)
#library(maps)
#library(colorRamps)
#library(scales)
#library(raster)
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
stats <- read.csv("FINALFIGS/2_basicStats/BasicStats_all.csv")
pca_results <- read.csv("FINALFIGS/PC.out+EnvData.csv")
#meta <- read.csv("../FINAL_DATASET/gigas_meta_41pop_edited.csv")

#pdf(width=9,height=4,"prettyHs.pdf")
#quartz(width=9,height=6)
par(mfrow=c(1,2),mar=c(4,5,1,1))

#### pc1 ~ sstmean
plot(PC1~BO_sstmean,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],xlim=range(pca_results$BO_sstmean),ylim=range(pca_results$PC1),col=pca_results$pc1.cols,pch=c(21,20)[factor(pca_results$natnon)],xlab="",ylab="")
abline(lm(PC1~BO_sstmean,data=pca_results))# & !rownames(pca_results)=="YOJ",]))
#pca_results$reg <- meta$GeneticRegions[match(pca_results2$Group.1,meta$pop)]
mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=pca_results$BO_sstmean,y=pca_results$PC1,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)#c("black","white")[factor(pca_results$natnon)])

pca_results$Hs <- stats$Hs[match(stats$X,pca_results$Group.1)]
plot(Hs~BO_sstmean,data=pca_results,cex=c(3,4)[factor(pca_results$natnon)],xlim=range(pca_results$BO_sstmean,na.rm=T),ylim=range(pca_results$Hs,na.rm=T),col=pca_results$pc1.cols,pch=c(21,20)[factor(pca_results$natnon)],xlab="",ylab="")
abline(lm(Hs~BO_sstmean,data=pca_results))# & !rownames(pca_results)=="YOJ",]))
#pca_results$reg <- meta$GeneticRegions[match(pca_results2$Group.1,meta$pop)]
mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Genetic PC1",side=2,line=2.5,cex=1.5)
text(x=pca_results$BO_sstmean,y=pca_results$PC1,pca_results$Group.1,cex=.5,col=pca_results$popTextCol)#c("black","white")[factor(pca_results$natnon)])





stats$lat <- pca_results$Latitude[match(stats$X,meta$pop)]
stats$lon <- meta$Longitude[match(stats$X,meta$pop)]
stats$sstmean <- extract(tmpmap,data.frame(stats$lon,stats$lat))
stats$sstmean[stats$X=="ALB"] <- 17.463# meta$satmean[meta$pop=="ALB"] <- 18.351
stats$sstmean[stats$X=="ESN"] <- 14.696#; meta$satmean[meta$pop=="ESN"] <- 12.403
# remove Chile (no lat lon) and MAI (n=4)
#stats <- stats[!stats$X%in%c("MAI"),]


stats$reg <- meta$Region[match(stats$X,meta$pop)]
stats$natnon <- ifelse(stats$reg%in%c("Japan","Korea"),"Native","Introduced")
#nat <- stats[natnon=="Native",]
#stats$GeneticRegions <- meta$GeneticRegions[match(stats$X,meta$pop)]
#stats$GeneticRegions[stats$natnon=="Introduced"] <- NA
#stats$GeneticRegions <- factor(stats$GeneticRegions)
#stats$GeneticRegions <- factor(stats$GeneticRegions,levels(stats$GeneticRegions)[c(1,4,6,5,2,3)])
#stats$GeneticRegions_col <- c("black", #1
#                             "red",#2
#                             "darkgreen",#3
#                             "darkgrey",#4
#                             "darkorchid",#5
#                             "deepskyblue")[factor(stats$GeneticRegions)]

#stats$GeneticRegions_col[stats$natnon=="Introduced"] <- "red"
stats$col <- pca_results$col[match(stats$X,pca_results$Group.1)]
stats <- stats[complete.cases(stats$Hs),]
plot(Hs~sstmean,data=stats,cex=c(3,4)[factor(stats$natnon)],ylim=range(stats$Hs),xlim=range(pca_results$sstmean),col=stats$col,pch=c(21,20)[factor(stats$natnon)],xlab="",ylab="")
abline(lm(Hs~sstmean,data=stats))# & !rownames(all)=="YOJ",]))
#all$reg <- meta$GeneticRegions[match(all2$Group.1,meta$pop)]
text(y=stats$Hs,x=stats$sstmean,stats$X,cex=.5,col=c("black","white")[factor(stats$natnon)])
mtext("Mean SST (ºC)",side=1,line=2.5,cex=1.5)
mtext("Hs",side=2,line=2.5,cex=1.5)
#mtext("PCA",line=-2)

print(cor.test(y=pca_results$PC1,x=pca_results$sstmean))
print(cor.test(y=stats$Hs,x=stats$sstmean))

stats$pc1 <- pca_results$PC1[match(stats$X,pca_results$Group.1)]
print(cor.test(y=stats$Hs,x=stats$pc1))

dev.off()

