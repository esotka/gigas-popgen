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
meta <- as.data.frame(read_xlsx("FINALFIGS/5_RandomForestCirclize/gverm/ece33001-sup-0002-tables1_EDITED.xlsx",sheet = 2))
meta <- meta[!meta$Subregion=="EUSA",]
## 1) find which 1º by 1º quadrat for each pop
tmp <- as.data.frame(meta)
gridIDs.tmp <- c()
for (i in 1:length(tmp$`Site abb.`))
{
  ppp.tmp <- ppp(x=tmp$Longitude[tmp$`Site abb.`==tmp$`Site abb.`[i]],y=tmp$Latitude[tmp$`Site abb.`==tmp$`Site abb.`[i]],xrange=domain[1:2],yrange=domain[3:4])
  grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
  gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$`Site abb.`)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$`Site abb.`,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
meta$sourceID <- ifelse(is.na(meta$sourceID),meta$Subregion,meta$sourceID)
tmpreg <- matrix(c(
  "bam","PNW",
  "bob","Cali",
  "eld","PNW",
  "elk","Cali",
  "ens","Cali",
  "ptw","PNW",
  "moo","PNW",
  "tmb","Cali"),nrow=8,ncol=2,byrow = T)
tmpreg2 <-tmpreg[match(meta$`Site abb.`,tmpreg[,1]),2]
meta$sourceID <- ifelse(is.na(tmpreg2),meta$sourceID,tmpreg2)
meta$country <- unlist(lapply(strsplit(meta$Site,", "),"[[",2))
meta$sourceID <- ifelse(meta$country%in%c("E","F","PT","MR"),"EuropeSouth",meta$sourceID)
meta$sourceID <- ifelse(meta$country%in%c("D","DK","IRE","UK"),"EuropeNorth",meta$sourceID)

### 2) generate the random forest - native pops vs introduced pops
dat <- read.csv('FINALFIGS/5_RandomForestCirclize/gverm/KruegerHadfieldetal_EcolEvol_diploidspsex_EDITED.csv',skip=2)[,1:22] # 10 usats
nloci = 10
colnames(dat) <- c("Ind","Pop",paste(sort(rep(letters[1:nloci],2)),rep(1:2,nloci),sep=""))
dat$Ind <- paste(1:nrow(dat),dat$Pop,sep="")
dat <- dat[dat$Pop%in%meta$`Site abb.`,] # remove east coast
sch = data.frame(pop=dat$Pop,source=meta$sourceID[match(dat$Pop,meta$`Site abb.`)])
sch$intro <- sch$source%in%c("EuropeNorth","EuropeSouth","PNW","Cali")

native_data = dat[!sch$intro,]
native_pops = as.factor(native_data$Pop)
intro_data =  dat[sch$intro,]
intro_pops =  as.factor(intro_data$Pop)

rf = randomForest(x=native_data,y=native_pops)

plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))

pdf("FINALFIGS/5_RandomForestCirclize/gverm/random_forest.pdf")
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()

write.csv(tbl,"FINALFIGS/5_RandomForestCirclize/gverm/RFprediction.csv",quote=F)

### 3) make circlize plot
rowReg <- meta$sourceID[match(rownames(tbl),meta$`Site abb.`)]
rowReg <- factor(rowReg); rowReg <- factor(rowReg,levels(rowReg)[c(1,2,6,5,3,4)])
colReg <- meta$sourceID[match(colnames(tbl),meta$`Site abb.`)]
colReg <- factor(colReg); colReg <- factor(colReg,levels(colReg)[c(1,4,2,3)])

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
#mat <- mat+.01

#cols.to.use <- c(meta2$pc1.cols[match(rownames(mat),meta2$GeneticRegions)],rep("grey",ncol(mat)))# rows first,  cols second
cols.to.use <- c(blue2red(5),"black",rep("grey",ncol(mat)))
#cols.to.use[is.na(cols.to.use)] <- "red"
rownames(mat) <- c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima","Korea / western Japan")


pdf("FINALFIGS/5_RandomForestCirclize/gverm/circlize.pdf",width=10,height=10)
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
mtext("G verm",line=-5,cex=2)
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
write.csv(mat,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/gvermByReg.csv")
