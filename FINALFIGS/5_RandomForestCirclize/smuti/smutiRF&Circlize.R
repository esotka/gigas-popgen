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
meta <- read.csv("FINALFIGS/5_RandomForestCirclize/smuti/LeCamTableS1_edited.csv")
## 1) find which 1º by 1º quadrat for each pop
tmp <- as.data.frame(meta)
gridIDs.tmp <- c()
for (i in 1:nrow(tmp))
{
  ppp.tmp <- ppp(x=tmp$Longitude[i],y=tmp$Latitude[i],xrange=domain[1:2],yrange=domain[3:4])
  if(ppp.tmp$n==0){gridIDs.tmp <- c(gridIDs.tmp,NA)}
  else{
    grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
    df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
    gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
  }
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$numpop)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]


meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$numpop,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
meta$sourceID <- ifelse(is.na(meta$sourceID),meta$Country,meta$sourceID)
meta <- meta[!meta$sourceID%in%c("Russia","Morocco","USA AK","Mexico","USA AK","USA OR"),]
tmpreg <- matrix(c(
  "USA,CA","Cali",
  "USA CA","Cali",
  "USA WA","PNW",
  "UK","EuropeNorth",
  "Spain","EuropeSouth",
  "Portugal","EuropeSouth",
  "Norway","EuropeNorth",
  "Ireland","EuropeNorth",
  "France","EuropeSouth",
  "Belgium","EuropeNorth"),nrow=10,ncol=2,byrow = T)
tmpreg2 <-tmpreg[match(meta$sourceID,tmpreg[,1]),2]
meta$sourceID <- ifelse(is.na(tmpreg2),meta$sourceID,tmpreg2)


### 2) generate the random forest - native pops vs introduced pops
dat <- read.delim('FINALFIGS/5_RandomForestCirclize/smuti/LeCam et al._Sargassum_Data&Microsatellite_edited.txt',skip=2) # 14 usats
nloci = 14
dat$pop <- gridIDs.meta$gridIDs[match(dat$pop,gridIDs.meta$pop)]
dat <- dat[!is.na(dat$pop),] # remove NAs ==> populations not in these grids
colnames(dat) <- c("Ind","Pop",paste(sort(rep(letters[1:nloci],2)),rep(1:2,nloci),sep=""))
dat <- dat[!is.na(meta$sourceID[match(dat$Pop,meta$gridIDs)]),]
dat <- dat[complete.cases(dat),]
sch = data.frame(pop=dat$Pop,source=meta$sourceID[match(dat$Pop,meta$gridIDs)])
sch$intro <- sch$source%in%c("EuropeNorth","EuropeSouth","PNW","Cali")

native_data = dat[!sch$intro,]
native_pops = as.factor(native_data$Pop)
intro_data =  dat[sch$intro,]
intro_pops =  as.factor(intro_data$Pop)

rf = randomForest(x=native_data,y=native_pops)

plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))

pdf("FINALFIGS/5_RandomForestCirclize/smuti/random_forest.pdf")
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
#heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_RandomForestCirclize/smuti/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$gridIDs)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels(rowReg)[c(1,3,2)])
colReg <- meta$sourceID[match(colnames(tbl),meta$gridIDs)]
colReg <- factor(colReg); colReg <- factor(colReg,levels(colReg)[c(1,4,2,3)])

dat2 <- as.matrix(tbl[order(rowReg),order(colReg)])

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
cols.to.use <- c(blue2red(5)[c(1,3,4)],rep("grey",ncol(mat)))
#cols.to.use[is.na(cols.to.use)] <- "red"
rownames(mat) <- c("Hokkaido",
                   #"Miyagi",
                   "Tokyo",
                   "Seto Inland Sea"
                   #"Kagoshima",
                   #"Korea / western Japan"
                   )


pdf("FINALFIGS/5_RandomForestCirclize/smuti/circlize.pdf",width=10,height=10)
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
mtext("Sarg muti",line=-5,cex=2)
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
write.csv(mat,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/smutiByReg.csv")
