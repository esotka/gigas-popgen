rm(list=ls())
library(spatstat)
library(readxl)
library(scales)
library(circlize)
library(colorRamps)
library(maps)
library(mapdata)
library(maptools)
library(ranger)
library(caret)
library(strataG)

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

### map ###
### map of populations ###
map1 <- function() {
  map("worldHires",xlim=c(117.5,150),ylim=c(30,47),col="gainsboro",fill=TRUE)
  tmp <- meta_source[!meta_source$sourceID=="nonSource",]
  rect(tmp$lon.sq1,tmp$lat.sq2,tmp$lon.sq2,tmp$lat.sq1,col=alpha(c("blue","darkgreen","lightgreen","red","black"),.5)[factor(tmp$sourceID)])
  text(x=c(144.7672,143.0341,140.8678,134.6165,122.67), y=c(41.37038, 38.15945, 34.36471,32,31),c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima"),col=c("blue","darkgreen","black","red","lightgreen"),pos=4,cex=.7)
  box()
  points(meta$Longitude,meta$Latitude,pch=20,cex=2)
}
pdf("FINALFIGS/5_RandomForestCirclize/gverm/map.pdf",height=6,width=5); map1(); dev.off()


### 2) generate the random forest - native pops vs introduced pops
dat <- read.csv('FINALFIGS/5_RandomForestCirclize/gverm/KruegerHadfieldetal_EcolEvol_diploidspsex_EDITED.csv',skip=2)[,1:22] # 10 usats
nloci = 10
colnames(dat) <- c("Ind","Pop",paste(sort(rep(letters[1:nloci],2)),rep(1:2,nloci),sep=""))
dat$Ind <- paste(1:nrow(dat),dat$Pop,sep="")
dat <- dat[dat$Pop%in%meta$`Site abb.`,] # remove east coast
sch = data.frame(pop=dat$Pop,source=meta$sourceID[match(dat$Pop,meta$`Site abb.`)])
sch$intro <- sch$source%in%c("EuropeNorth","EuropeSouth","PNW","Cali")

gi <- df2gtypes(dat, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
all_mat = as.data.frame(gtypes2genind(gi)@tab)
native_data = all_mat[!sch$intro,]
native_pops = sch$pop[!sch$intro]
native_reg = as.factor(meta$sourceID[match(native_pops,meta$`Site abb.`)])
native_reg = factor(native_reg,levels=levels(native_reg)[c(1,2,6,5,3,4)])
intro_data = all_mat[sch$intro,]
intro_pops = sch$pop[sch$intro]
intro_reg = meta$sourceID[match(intro_pops,meta$`Site abb.`)]
intro_reg <- as.factor(intro_reg)

native = data.frame(cbind(native_reg,native_data))
names(native)[1]="reg"
intro = data.frame(cbind(intro_pops,intro_data))
names(intro)[1]="pop"


### RF ####
########################## training and testing

save(native,file="FINALFIGS/5_RandomForestCirclize/crossValidation/GvermforRangerCrossValidation.rda")

#fitControl=trainControl(method="repeatedcv", number=10,repeats=10) #10-fold cv repeated 10 times
#rangerGrid = expand.grid(mtry = round(seq(1000,2000,length=5)),splitrule=c("gini","extratrees"),min.node.size=c(1,3,5,10) )
#rangerGrid = expand.grid(mtry = round(seq(2,20,length=10)),splitrule=c("gini","extratrees"),min.node.size=c(1,3,5,10) )

rangerGrid = expand.grid(mtry = 150,splitrule=c("extratrees"),min.node.size=c(1) )

#control = list(ranger=list(method="ranger",tuneGrid=rangerGrid))

#rfits = lapply(control,function(x)
#{
#  print(x$method)
#  train(reg~., data=train, method=x$method,
#        tuneGrid=x$tuneGrid, trControl=fitControl)
#})

#rf = ranger(pop~.,data=native, mtry=rfits[[1]]$bestTune$mtry,
#           splitrule=rfits[[1]]$bestTune$splitrule,
#            min.node.size=rfits[[1]]$bestTune$min.node.size)

rf = ranger(reg~.,data=native, mtry=rangerGrid[[1]],
            splitrule=rangerGrid[[2]],
            min.node.size=rangerGrid[[3]])


pred = predict(rf,data=intro)$predictions
actual=intro_reg

tbl = table(pred,actual)
save(file="FINALFIGS/5_RandomForestCirclize/gverm/region_assignment.rda",tbl)
pdf("FINALFIGS/5_RandomForestCirclize/gverm/region_assignment.pdf")
par(mar=c(6,4,2,2))
heatmap(tbl)
dev.off()

### circlize plot
mat <- as.matrix(tbl)
mat <- mat+.01

cols.to.use <- c(blue2red(5),"black",rep("grey",ncol(mat)))
#rownames(mat) <- c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima","Korea / western Japan")


pdf("FINALFIGS/5_RandomForestCirclize/gverm/region_assignment_circlize.pdf",width=10,height=10)
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

## write sample sizes for summary
n <- data.frame(n=c(table(native_pops),table(intro_pops)))
n$reg <- meta$sourceID[match(rownames(n),meta$pop)]
write.csv(n,"FINALFIGS/5_RandomForestCirclize/ALLSPECIES/cgigas_sampleSize.csv",quote=F)
