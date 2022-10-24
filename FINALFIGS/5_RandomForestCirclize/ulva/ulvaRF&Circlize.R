rm(list=ls())
library(randomForest)
library(readxl,warn.conflicts = F,quietly = T)
library(spatstat,warn.conflicts = F,quietly = T) # marks, ppp, quadratcount
library(scales,warn.conflicts = F,quietly = T) # alpha
library(circlize)
library(colorRamps)
library(reshape)

domain <- c(-130,180,-50,75) 
load("FINALFIGS/0_globalGrid/df.globe.Rda")
meta_source <- read.csv("FINALFIGS/0_globalGrid/df.globe_source.csv")
meta <- read.delim("FINALFIGS/5_RandomForestCirclize/ulva/Hanyuda et al 2018 Supp TableS6 Ulva.txt")
meta <- unique(meta[,c("Name1","lat","lon","Country")])
## 1) find which 1º by 1º quadrat for each pop
tmp <- as.data.frame(meta)
gridIDs.tmp <- c()
for (i in 1:length(tmp$Name1))
{
  ppp.tmp <- ppp(x=tmp$lon[i],y=tmp$lat[i],xrange=domain[1:2],yrange=domain[3:4])
  if(ppp.tmp$n==0){gridIDs.tmp <- c(gridIDs.tmp,NA)}
  else{
    grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
    df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
    gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
  }}

gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$Name1)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$Name1,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
meta$sourceID <- ifelse(is.na(meta$sourceID),meta$Country,meta$sourceID)

## remove Japan pop that wasn't on a coast
meta <- meta[!meta$sourceID%in%c("Japan","Korea","Austalia","Australia","Mexico"),]
meta$sourceID[meta$sourceID%in%c("France")] <- "EuropeSouth"
meta$sourceID[meta$sourceID%in%c("Netherlands")] <- "EuropeNorth"

### 2) generate the random forest - native pops vs introduced pops
dat <- read.delim("FINALFIGS/5_RandomForestCirclize/ulva/Hanyuda et al 2018 Supp TableS6 Ulva.txt")[,c("Name1","Haplotype","Number")]
#dat$sourceID <- meta$sourceID[match(dat$Name1,meta$Name1)]
dat <- dat[dat$Name1%in%meta$Name1,]
#md <- melt(dat,id="Haplotype")
md <- melt(dat,id=c("Name1","Haplotype"))
#md2 <- cast(md,Haplotype~Locality.codea) # hap = row; pop=col
#rownames(md2) <- md2$Haplotype
#md2 <- md2[,-1]

#md <- md[complete.cases(md),]
popID <- rep(as.character(md$Name1),md$value)
hapInd <- rep(md$Haplotype,md$value)
indID <- paste(popID,1:length(hapInd),sep="_")
metaInd <- data.frame(popID,hapInd)
md2 <- as.matrix(table(metaInd))
md2_pop <- rownames(md2)#unlist(lapply(strsplit(rownames(md2),"_"),"[[",1))
md2_source <- meta$sourceID[match(md2_pop,meta$Name1)]

native_data = md2[!md2_source%in%c("NZ","USA","EuropeNorth","EuropeSouth"),]
native_pops = as.factor(md2_pop[!md2_source%in%c("NZ","USA","EuropeNorth","EuropeSouth")])
intro_data =  md2[md2_source%in%c("NZ","USA","EuropeNorth","EuropeSouth"),]
intro_pops =  as.factor(md2_pop[md2_source%in%c("NZ","USA","EuropeNorth","EuropeSouth")])


rf = randomForest(x=native_data,y=native_pops)

plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))

pdf("FINALFIGS/5_RandomForestCirclize/ulva/random_forest.pdf")
#par(mar=c(7,4,4,4)+0.1)
heatmap(tbl)#,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,])#,cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_RandomForestCirclize/ulva/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$Name1)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels=levels(rowReg)[c(1,2,6,5,3,4)])
colReg <- meta$sourceID[match(colnames(tbl),meta$Name1)]
colReg[colReg=="USA"] <- "NAm_south"
colReg <- factor(colReg)#; colReg <- factor(colReg,levels(colReg)[c()])

#tbl <- tbl[order(rowReg),order(colReg)]
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

rownames(datByReg2) <- c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima","Korea / wJapan")
mat <- as.matrix(datByReg2)
mat <- mat+0.01
cols.to.use <- c(blue2red(5),"black",rep("grey",ncol(mat)))

pdf("FINALFIGS/5_RandomForestCirclize/ulva/circlize.pdf",width=10,height=10)
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
mtext("Ulva",line=-5,cex=2)
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
write.csv(mat,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/ulvaByReg.csv")

## write sample sizes for summary
n <- data.frame(n=c(table(native_pops),table(intro_pops)))
n$reg <- c(as.character(rowReg),as.character(colReg))
write.csv(n,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/ulva_sampleSize.csv",quote=F)
