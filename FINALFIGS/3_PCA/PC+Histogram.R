### points of PC1 output of SNPs
rm(list=ls())
library(maps);library(fields);library(reshape);library(mapdata);library(colorRamps);library(plotrix);library(scales);library(ggplot2);library(ggridges); library(gridExtra)

### PCA ####
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
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

### PCA on everything
xbar <- aggregate(pca1_mat[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$Region2[match(xbar$Group.1,meta$pop)]
cols.to.use <- c(blue2red(5),"black")
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima","nonSource")
xbar$pc1.cols = cols.to.use[match(xbar$GeneticRegions,names(cols.to.use))]
xbar$popTextCol = "black"
xbar$popTextCol[xbar$Group.1%in%c("SAR","GOS","GWG","KOJ","MAI","AKK")] <- "white"

pops = substr(rownames(pca1_mat),1,3)
pts.cols = xbar$pc1.cols[match(pops,xbar$Group.1)]
segments_df = data.frame(pca1_mat,xbar_pc1 = xbar$PC1[match(pops,xbar$Group.1)],xbar_pc2 = xbar$PC2[match(pops,xbar$Group.1)])

f1 = ggplot() +
     geom_point(data=pca1_mat,aes(x = PC2, y = PC1),color = alpha(pts.cols,0.2)) + 
     geom_point(data = xbar,aes(x=PC2,y=PC1),color= xbar$pc1.cols,size=6) + 
     geom_segment(data=segments_df,aes(x=PC2,y=PC1,xend=xbar_pc2,yend=xbar_pc1),color=alpha(pts.cols,0.2)) +
     #scale_y_continuous(limits=c(-10,10)) + 
     theme_light() +
     annotate(geom = "text", x = xbar$PC2, y = xbar$PC1,size=2, label = xbar$Group.1,col=xbar$popTextCol)

# histogram
# remove YOJ
pca1_mat_sub = pca1_mat[!nat.pop=="YOJ",]
nat.pop_sub = factor(substr(rownames(pca1_mat_sub),1,3))
nat.pop_sub = factor(nat.pop_sub,levels(nat.pop_sub)[c(15,1,8,9,10,14,2,11,5,12,13,16,3,4,6,7)])
f4 = ggplot(pca1_mat_sub,aes(x=PC1,y=factor(nat.pop_sub),fill=factor(nat.pop_sub))) + 
     geom_density_ridges(scale=1,alpha=2) + 
     scale_fill_cyclical(values = xbar$pc1.cols[match(levels(nat.pop_sub),xbar$Group.1)])+
     theme_minimal()#+ 



### PC and assignment on Akkeshi (SAR only) + Miyagi, Seto and Tokyo (remove AKK and YOJ) ###
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
pop_nat = meta$pop[meta$Region2 %in% c("Tokyo","Seto","Miyagi") | meta$pop=="SAR"]
geno <- read.delim("FINALFIGS/SNPs_noZeros.txt",header=T)
rownames(geno) <- geno$id
geno <- geno[,-1]
pop <- substr(rownames(geno),1,3)
nat <- geno[pop%in%pop_nat,]
nat.pop <- factor(substr(rownames(nat),1,3))
nat.reg <- factor(meta$Region2[match(nat.pop,meta$pop)])
nat.reg <- factor(nat.reg,levels(nat.reg)[c(1,2,4,3)])

pca1 <- prcomp(nat)
pca1_mat = pca1$x[,1:2]

### PCA # 2
xbar <- aggregate(pca1_mat[,1:2],by=list(nat.pop),mean)
xbar$GeneticRegions <- meta$Region2[match(xbar$Group.1,meta$pop)]
cols.to.use <- c(blue2red(5)[-5])
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto")
xbar$pc1.cols = cols.to.use[match(xbar$GeneticRegions,names(cols.to.use))]
xbar$popTextCol = "black"
xbar$popTextCol[xbar$Group.1%in%c("SAR","GOS","GWG","KOJ","MAI","AKK")] <- "white"

pops = substr(rownames(pca1_mat),1,3)
pts.cols = xbar$pc1.cols[match(pops,xbar$Group.1)]
segments_df = data.frame(pca1_mat,xbar_pc1 = xbar$PC1[match(pops,xbar$Group.1)],xbar_pc2 = xbar$PC2[match(pops,xbar$Group.1)])

f2 = ggplot() +
     geom_point(data=pca1_mat,aes(x = PC2, y = -PC1),color = alpha(pts.cols,0.2)) + 
     geom_point(data = xbar,aes(x=PC2,y=-PC1),color= xbar$pc1.cols,size=6) + 
     geom_segment(data=segments_df,aes(x=PC2,y=-PC1,xend=xbar_pc2,yend=-xbar_pc1),color=alpha(pts.cols,0.2)) +
     theme_light() +
     annotate(geom = "text", x = xbar$PC2, y = -xbar$PC1,size=2, label = xbar$Group.1,col=xbar$popTextCol)

pop_non <- meta$pop[meta$NatNon=="Introduced"]
non = geno[pop%in%pop_non,]
non.pop <- substr(rownames(non),1,3)
non.reg <- meta$Region2[match(non.pop,meta$pop)]

introPredict4 <- predict(pca1,newdata=non)

df = data.frame(pop=c(as.character(nat.pop),non.pop),reg=c(as.character(nat.reg),non.reg),pc1=-c(pca1$x[,1],introPredict4[,1]))
## remove Kag and separate AKK and SAR in PC
#natnon <- natnon[!(natnon$reg%in%c("nonSource","Kagoshima","Akkeshi","Miyagi","Seto","Tokyo")),]
df$reg[df$reg == "Akkeshi"] = "Hokkaido"
df$reg[df$reg=="Seto"] = "Seto Inland Sea"
df$reg[df$reg=="soCalifornia"] = "N America (Cali)"
df$reg[df$reg=="PNW"] = "N America (CA/WA)"
df$reg[df$reg=="noEurope"] = "N Europe (UK/De/G/Sw)"
df$reg[df$reg=="soEurope"] = "S Europe (Fr/Sp)"

df$reg <- factor(df$reg)
df$reg <- factor(df$reg,levels=levels(df$reg)[c(3,4,11,10,8,6,5,1,2,9,7)])
f3= ggplot(df,aes(x = pc1, y = reg, fill = reg)) + 
     geom_density_ridges(scale=1,alpha=2) + 
     scale_fill_cyclical(values = as.character(c(cols.to.use,rep("grey",7)))) +
     theme_minimal()

#quartz()
png("FINALFIGS/3_PCA/PC-all.png",width=5,height=5,units="in",res=400); f1; dev.off()

png("FINALFIGS/3_PCA/PC-histograms_japan.png",width=5,height=5,units="in",res=400); f4; dev.off()

png("FINALFIGS/3_PCA/PC+Histogram.png",width=9,height=5,units="in",res=400)
grid.arrange(grobs = list(f2,f3),nrow=1,padding=0)#layout_matrix = c(1,2,2))
dev.off()

