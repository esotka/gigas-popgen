### points of PC1 output of SNPs
rm(list=ls())
library(maps);library(fields);library(reshape);library(mapdata);library(colorRamps);library(plotrix);library(scales);library(ggplot2);library(ggridges); library(gridExtra)
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
pca1_mat = pca1$x[,1:2]
### predict ###
introPredict4 <- predict(pca1,newdata=non)

natnon <- data.frame(pop=c(as.character(nat.pop),non.pop),reg=c(as.character(nat.reg),non.reg),pc1=c(pca1$x[,1],introPredict4[,1]))
## remove Kag and non-source
natnon <- natnon[!(natnon$reg%in%c("nonSource","Kagoshima","Akkeshi","Miyagi","Seto","Tokyo")),]
natnon$reg <- factor(natnon$reg)
natnon$reg <- factor(natnon$reg,levels=levels(natnon$reg)[c(3,6,5,1,2,4,7)])
f1= ggplot(natnon,aes(x = pc1, y = reg, fill = stat(x))) + 
    geom_density_ridges_gradient(scale=1,show.legend = F,rel_min_height = 0.01) + 
    scale_x_continuous(limits=c(-10,10)) + 
    theme_minimal() + 
    coord_flip() 
   #scale_fill_viridis_c(name = "PC1", option = "C") + 
 f1
#quartz(width=9,height=5)
#par(mar=c(1,1,1,1),mfrow=c(2,1))


# layout
#png("FINALFIGS/3_PCA/PC+Histogram.png",width=6,height=6,units="in",res=500); plot(f1); dev.off()

#png("FINALFIGS/3_PCA/PC_globe_separateNatNon.png",width=6,height=6,units="in",res=500)
#par(mar=c(2,4,2,4))
xbar <- aggregate(pca1_mat[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$Region2[match(xbar$Group.1,meta$pop)]
cols.to.use <- c(blue2red(5),"black")
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima","nonSource")
xbar$pc1.cols = cols.to.use[match(xbar$GeneticRegions,names(cols.to.use))]
xbar$popTextCol = "black"
xbar$popTextCol[xbar$Group.1%in%c("SAR","GOS","GWG","KOJ","MAI","AKK")] <- "white"

pops = substr(rownames(pca1_mat),1,3)
pts.cols = xbar$pc1.cols[match(pops,xbar$Group.1)]
#popPts <- xbar[match(pops,xbar$Group.1),]
segments_df = data.frame(pca1_mat,xbar_pc1 = xbar$PC1[match(pops,xbar$Group.1)],xbar_pc2 = xbar$PC2[match(pops,xbar$Group.1)])
#seg_col <- xbar$pc1.cols[match(nat.pop,as.character(xbar$Group.1))]
#segments(x0 = -pca1$x[,2],y0=pca1$x[,1],x1=-popPts$PC2,y1=popPts$PC1,col=alpha(seg_col,0.2),lwd=2)#,col=alpha(col.sub,.5))


f2 = ggplot() +
     geom_point(data=pca1_mat,aes(x = PC2, y = PC1),color = alpha(pts.cols,0.2)) + 
     geom_point(data = xbar,aes(x=PC2,y=PC1),color= xbar$pc1.cols,size=6) + 
     geom_segment(data=segments_df,aes(x=PC2,y=PC1,xend=xbar_pc2,yend=xbar_pc1),color=alpha(pts.cols,0.2)) +
     scale_y_continuous(limits=c(-10,10)) + 
     theme_light() +
     annotate(geom = "text", x = xbar$PC2, y = xbar$PC1,size=2, label = xbar$Group.1,col=xbar$popTextCol)
#f2
#quartz(width=9,height=5)

png("FINALFIGS/3_PCA/PC+Histogram.png",width=9,height=5,units="in",res=400)
grid.arrange(f2,f1,ncol=2,nrow=1)#layout = c(1,2))

dev.off()

