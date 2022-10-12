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
pdf("FINALFIGS/5_RandomForestCirclize/cgiga/random_forest.pdf")
plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_RandomForestCirclize/cgiga/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$pop)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels(rowReg)[c(1,2,6,5,3,4)])
colReg <- meta$sourceID[match(colnames(tbl),meta$pop)]
colReg[colReg%in%c("France","Spain","Ireland","Sweden","Norway","Denmark")] <- "Europe"
colReg <- factor(colReg); colReg <- factor(colReg,levels(colReg)[c(4,1:2,6,5,3)])

dat2 <- as.matrix(tbl[order(rowReg),order(colReg)])
df = data.frame(from = rep(rownames(dat2), times = ncol(dat2)),
                to = rep(colnames(dat2), each = nrow(dat2)),
                value = as.vector(dat2),
                stringsAsFactors = FALSE)
## by reg


datByReg <- c()
for (i in 1:length(levels(rowReg))) # native
{
  tmp <- tbl[rowReg==levels(rowReg)[i],]
  if(is.null(nrow(tmp)))
  {datByReg <- rbind(datByReg,tmp)}
  else{
    datByReg <- rbind(datByReg,colSums(tbl[rowReg==levels(rowReg)[i],]))
  }}
datByReg2 <- c()
for(j in 1:length(levels(colReg)))
{  
  tmp <- datByReg[,colReg==levels(colReg)[j]]
  if(is.null(ncol(tmp)))
  {datByReg2 <- cbind(datByReg2,tmp)}
  else{
    datByReg2 <- cbind(datByReg2,rowSums(tmp))
  }}

rownames(datByReg2) <- levels(rowReg)
colnames(datByReg2) <- levels(colReg)
#datByReg2 <- datByReg2[,c(2:7,1)]
#datByReg2 <- rbind(datByReg2,rep(0,ncol(datByReg2))) # add dummy region
mat <- as.matrix(datByReg2)
mat <- mat+.01

#cols.to.use <- c(meta2$pc1.cols[match(rownames(mat),meta2$GeneticRegions)],rep("grey",ncol(mat)))# rows first,  cols second
cols.to.use <- c(blue2red(5),"black",rep("grey",ncol(mat)))
#cols.to.use[is.na(cols.to.use)] <- "red"
rownames(mat) <- c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima","Korea / western Japan")


pdf("FINALFIGS/5_RandomForestCirclize/cgiga/circlize.pdf",width=10,height=10)
circos.clear()
circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)

chordDiagram(x = mat, 
             directional = 1,
             grid.col = cols.to.use,
             annotationTrack = "grid", 
             transparency = 0.25,  
             annotationTrackHeight = c(0.1, 0.1),
             #diffHeight  = -0.04,
             direction.type = c("arrows", "diffHeight"),
             link.arr.type = "big.arrow",
             link.arr.length =  0.15,
             diffHeight = -0.001,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
mtext("C gigas",line=-5,cex=2)
# we go back to the first track and customize sector labels
#circos.track(track.index = 1, panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
#}, bg.border = NA) # here set bg.border to NA is important
circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  if(abs(xplot[2] - xplot[1]) < 10) {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",niceFacing = TRUE, adj = c(0, 0.5), col = "black",cex=1.5)
  } else {
    circos.text(mean(xlim), ylim[1], sector.name, facing = "inside", niceFacing = TRUE, adj = c(0.5, 0), col= "black",cex=1.5)
  }
}, bg.border = NA)

dev.off()
write.csv(mat,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/cgigasByReg.csv")
