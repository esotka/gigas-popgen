rm(list=ls())
library(readxl,warn.conflicts = F,quietly = T)
library(spatstat,warn.conflicts = F,quietly = T) # marks, ppp, quadratcount
library(scales,warn.conflicts = F,quietly = T) # alpha
library(circlize)
library(colorRamps)
library(reshape)
library(maps)
library(mapdata)
library(ranger)


domain <- c(-130,180,-50,75) 
load("../df.globe.Rda")
meta_source <- read.csv("../df.globe_source.csv")
meta <- read.csv("Uwai_211inds.meta.csv")
meta <- unique(meta[,c("popID","lonInd","latInd","country")])
## 1) find which 1º by 1º quadrat for each pop
tmp <- as.data.frame(meta)
gridIDs.tmp <- c()
for (i in 1:length(tmp$popID))
{
  ppp.tmp <- ppp(x=tmp$lonInd[tmp$popID==tmp$popID[i]],y=tmp$latInd[tmp$popID==tmp$popID[i]],xrange=domain[1:2],yrange=domain[3:4])
  grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
  gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$popID)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

meta$gridIDs <- gridIDs.meta$gridIDs[match(meta$popID,gridIDs.meta$pop)]
meta$sourceID <- meta_source$sourceID[match(meta$gridIDs,meta_source$gridID)]
meta$sourceID <- ifelse(is.na(meta$sourceID),meta$country,meta$sourceID)

### map of populations ###
map1 <- function() {
  map("worldHires",xlim=c(117.5,150),ylim=c(30,47),col="gainsboro",fill=TRUE)
  tmp <- meta_source[!meta_source$sourceID=="nonSource",]
  tmp$sourceID <- factor(tmp$sourceID)
  tmp$sourceID <- factor(tmp$sourceID,levels=levels(tmp$sourceID)[c(1,2,5,4,3)])
  rect(tmp$lon.sq1,tmp$lat.sq2,tmp$lon.sq2,tmp$lat.sq1,col=alpha(blue2red(5),.5)[tmp$sourceID])
  text(x=c(144.7672,143.0341,140.8678,134.6165,122.67), y=c(41.37038, 38.15945, 34.36471,32,31),c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima"),col=blue2red(5),pos=4,cex=.7)
  box()
  points(meta$lonInd,meta$latInd,pch=20,cex=2)
}
pdf("map.pdf",height=6,width=5); map1(); dev.off()
png("map.png",height=5,width=5,units="in",res=400); map1(); dev.off()


### 2) generate the random forest - native pops vs introduced pops
dat <- read.csv("Uwai_211inds.meta.csv")
metaInd <- data.frame(indID=dat$indID,hapInd=dat$hapInd)
md2 <- as.matrix(table(metaInd))
md2_pop <- dat$popID[match(rownames(md2),dat$indID)]#unlist(lapply(strsplit(rownames(md2),"_"),"[[",1))
md2_source <- meta$sourceID[match(md2_pop,meta$popID)]

native_data = md2[!md2_source%in%c("New Zealand","France","USA"),]
native_pops = as.factor(md2_pop[!md2_source%in%c("New Zealand","France","USA")])
native_reg = md2_source[!md2_source%in%c("New Zealand","France","USA")]
intro_data =  md2[md2_source%in%c("New Zealand","France","USA"),]
intro_pops =  as.factor(md2_pop[md2_source%in%c("New Zealand","France","USA")])
intro_reg = md2_source[md2_source%in%c("New Zealand","France","USA")]

native = data.frame(reg=native_reg,cbind(native_data))
native = native[!as.character(native$reg)%in%c("China","Japan"),]
native$reg = factor(native$reg)
native$reg = factor(native$reg,levels=levels(native$reg)[c(1,2,6,5,3,4)])
intro = data.frame(pop=intro_pops,cbind(intro_data))


### RF ####
########################## training and testing

save(native,file="upinnForRangerCrossValidation.rda")

### ran on "XXXXCrossvalidation_forServer.R" separately 
load(file="upinn_Ranger_out.rda")

rf = ranger(reg~.,data=native, mtry=rfits[[1]]$bestTune$mtry,splitrule=rfits[[1]]$bestTune$splitrule,min.node.size=rfits[[1]]$bestTune$min.node.size)
print(rfits[[1]]$bestTune)
print(accuracy <- rfits[[1]]$results$Accuracy[as.numeric(rownames(rfits[[1]]$bestTune))])

pred = predict(rf,data=intro)$predictions
actual=intro_reg

tbl = table(pred,actual)
save(file="region_assignment.rda",tbl)
pdf("region_assignment.pdf")
par(mar=c(6,4,2,2))
heatmap(tbl)
dev.off()

write.csv(tbl,"RFprediction.csv",quote=F)

### 3) make circlize plot

mat <- as.matrix(tbl)
mat <- mat+0.01
mat = mat[,-2] # no NZ
cols.to.use <- c(blue2red(5),"black","black","darkgrey")

pdf("circlize_noLabels.pdf",width=5,height=5)
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
dev.off()
pdf("circlize.pdf",width=10,height=10)
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
write.csv(mat,"upinnByReg.csv")

## write sample sizes for summary
n <- data.frame(n=c(table(native_pops),table(intro_pops)))
n$reg <- meta$sourceID[match(rownames(n),meta$popID)]
n = n[!n$reg%in%c("China","Japan"),]
write.csv(n,"upinn_sampleSize.csv",quote=F)
