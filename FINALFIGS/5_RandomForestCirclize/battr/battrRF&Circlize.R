rm(list=ls())
library(randomForest)
library(readxl,warn.conflicts = F,quietly = T)
library(spatstat,warn.conflicts = F,quietly = T) # marks, ppp, quadratcount
library(scales,warn.conflicts = F,quietly = T) # alpha
library(circlize)
library(colorRamps)

domain <- c(-130,180,-50,75) 
load("FINALFIGS/0_globalGrid/df.globe.Rda")
meta_source <- read.csv("FINALFIGS/0_globalGrid/df.globe_source.csv")
meta <- as.data.frame(read_xlsx("FINALFIGS/5_RandomForestCirclize/battr/Miura2006haplotypes_meta.xlsx",sheet = 1))
## 1) find which 1º by 1º quadrat for each pop
tmp <- as.data.frame(meta)
gridIDs.tmp <- c()
for (i in 1:length(tmp$Site))
{
  ppp.tmp <- ppp(x=tmp$longitude[tmp$Site==tmp$Site[i]],y=tmp$latitude[tmp$Site==tmp$Site[i]],xrange=domain[1:2],yrange=domain[3:4])
  grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
  gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$Site)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$Site,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
tmpreg <- matrix(c(
  "Boundary","PNW",
  "Padila","PNW",
  "Bolinus","Cali",
  "Elkhorn","Cali"),nrow=4,ncol=2,byrow = T)
tmpreg2 <-tmpreg[match(meta$Site,tmpreg[,1]),2]
meta$sourceID <- ifelse(is.na(tmpreg2),meta$sourceID,tmpreg2)


### 2) generate the random forest - native pops vs introduced pops
dat <- as.data.frame(read_xlsx('FINALFIGS/5_RandomForestCirclize/battr/Miura2006haplotypes_haps.xlsx',sheet=1))
rownames(dat) <- dat$Haplotype; dat <- dat[,-1]
mat <- as.matrix(dat)
mat[is.na(mat)] <- 0
mat <- t(mat) ### row = pop; col = haplo

native_data = mat[1:14,]
native_pops = as.factor(rownames(mat)[1:14])
intro_data =  mat[15:18,]
intro_pops =  as.factor(rownames(mat)[15:18])


rf = randomForest(x=native_data,y=native_pops)

plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))

pdf("FINALFIGS/5_RandomForestCirclize/battr/random_forest.pdf")
#par(mar=c(7,4,4,4)+0.1)
heatmap(tbl)#,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,])#,cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_RandomForestCirclize/battr/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$Site)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels(rowReg)[c(1,5,4,2,3)])
colReg <- meta$sourceID[match(colnames(tbl),meta$Site)]
colReg <- factor(colReg)#; colReg <- factor(colReg,levels(colReg)[c()])

#dat <- as.data.frame(tbl[order(rowReg),order(colReg)])
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

rownames(datByReg2) <- c("Miyagi","Tokyo","Seto Inland Sea","Kagoshima","Korea / westernJapan")
mat <- as.matrix(datByReg2)
mat <- mat+0.01
cols.to.use <- c(blue2red(5)[2],rep("white",3),"white",rep("grey",ncol(mat)))

pdf("FINALFIGS/5_RandomForestCirclize/battr/circlize.pdf",width=10,height=10)
circos.clear()
circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)

chordDiagram(x = mat, 
             directional = 1,
             grid.col = cols.to.use,
             annotationTrack = "grid", 
             transparency = 0.25,  
             annotationTrackHeight = c(0.1, 0.1),
             #diffHeight  = -0.04,
             #keep.diagonal=T,
             #reduce = -1,
             self.link=10,
             direction.type = c("arrows", "diffHeight"),
             link.arr.type = "big.arrow",
             link.arr.length =  0.15,
             diffHeight = -0.001,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
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
