rm(list=ls())
library(randomForest)
library(spatstat,warn.conflicts = F,quietly = T) # marks, ppp, quadratcount
library(scales,warn.conflicts = F,quietly = T) # alpha
library(circlize)
library(RColorBrewer)
library(colorRamps)

domain <- c(-130,180,-50,75) 
load("FINALFIGS/0_globalGrid/df.globe.Rda")
meta_source <- read.csv("FINALFIGS/0_globalGrid/df.globe_source.csv")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")

## 1) find which 1º by 1º quadrat for each pop
tmp <- meta
gridIDs.tmp <- c()
for (i in 1:length(tmp$pop))
{
  ppp.tmp <- ppp(x=tmp$Longitude[tmp$pop==tmp$pop[i]],y=tmp$Latitude[tmp$pop==tmp$pop[i]],xrange=domain[1:2],yrange=domain[3:4])
  grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
  gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$pop)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$pop,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
meta$sourceID <- ifelse(is.na(meta$sourceID),meta$Region,meta$sourceID)



### 2) generate the random forest - native pops vs introduced pops
tab = read.table("FINALFIGS/SNPs_noZeros.txt",header=T)
rownames(tab)=tab[,1]
tab=tab[,-1]

pc = prcomp(tab)
tab.pred=predict(pc)
tab.pred = tab
pop=sapply(strsplit(rownames(tab),"_"),function(x)x[1])

native_data = tab.pred[pop %in% meta$pop[meta$NatNon=="Native"],]
native_pops = as.factor(pop[pop %in% meta$pop[meta$NatNon=="Native"]])
intro_data =  tab.pred[!(pop %in% meta$pop[meta$NatNon=="Native"]),]
intro_pops =  as.factor(pop[!(pop %in% meta$pop[meta$NatNon=="Native"])])

rf = randomForest(x=native_data,y=native_pops)
pdf("FINALFIGS/5_circlize/cgiga/random_forest.pdf")
plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_circlize/cgiga/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$pop)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels(rowReg)[c(1,2,6,5,3,4)])
colReg <- meta$sourceID[match(colnames(tbl),meta$pop)]
colReg <- factor(colReg); colReg <- factor(colReg,levels(colReg)[c(1:2,6,9,8,10,4,3,5,7,11)])

dat2 <- as.matrix(tbl[order(rowReg),order(colReg)])
df = data.frame(from = rep(rownames(dat2), times = ncol(dat2)),
                to = rep(colnames(dat2), each = nrow(dat2)),
                value = as.vector(dat2),
                stringsAsFactors = FALSE)
