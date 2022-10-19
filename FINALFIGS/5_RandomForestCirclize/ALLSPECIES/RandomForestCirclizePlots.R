# make a PDF with all chordDiagrams
# make sure you run all of the species-specific RF files first (e.g., rudiRF&Circlize.R)
library(circlize)
library(colorRamps)
rm(list=ls())
giga <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/cgigasByReg.csv"); rownames(giga) <- giga[,1]; giga <- as.matrix(giga[,-1]) 
gverm <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/gvermByReg.csv"); rownames(gverm) <- gverm[,1]; gverm <- as.matrix(gverm[,-1])
batt <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/battrByReg.csv"); rownames(batt) <- batt[,1]; batt <- as.matrix(batt[,-1])
batt_HL1 <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/battr_HL1_ByReg.csv"); rownames(batt_HL1) <- batt_HL1[,1]; batt_HL1 <- as.matrix(batt_HL1[,-1])
batt_HL6 <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/battr_HL6_ByReg.csv"); rownames(batt_HL6) <- batt_HL6[,1]; batt_HL6 <- as.matrix(batt_HL6[,-1])
upinn <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/upinnByReg.csv"); rownames(upinn) <- upinn[,1]; upinn <- as.matrix(upinn[,-1]); # remove Mexico
smuti <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/smutiByReg.csv"); rownames(smuti) <- smuti[,1]; smuti <- as.matrix(smuti[,-1])
cutl.tmp <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/cutlByReg.csv"); cutl <- as.matrix(cutl.tmp[,-1]); rownames(cutl) <- cutl.tmp[,1]; colnames(cutl) <- colnames(cutl.tmp)[-1]
cera.tmp <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/ceraByReg.csv"); cera <- as.matrix(cera.tmp[,-1]); rownames(cera) <- cera.tmp[,1]; colnames(cera) <- colnames(cera.tmp)[-1]
rudi <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/rudiByReg.csv"); rownames(rudi) <- rudi[,1]; rudi <- as.matrix(rudi[,-1])
ulva <- read.csv("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/ulvaByReg.csv"); rownames(ulva) <- ulva[,1]; ulva <- as.matrix(ulva[,-1])

spp <- list(giga=giga,
            gverm=gverm,
            batt=batt,
            batt_HL1=batt_HL1,
            batt_HL6=batt_HL6,
            upinn=upinn,
            smuti=smuti,
            cutl=cutl,
            cera=cera,
            rudi=rudi,
            ulva=ulva)
cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(spp[["giga"]]))),
  gverm=c(blue2red(5),"black",rep("grey",ncol(spp[["gverm"]]))),
  batt=c(blue2red(5)[-1],"black",rep("grey",ncol(spp[["batt"]]))),
  batt_HL1=c(blue2red(5)[-1],"black",rep("grey",ncol(spp[["batt_HL1"]]))),
  batt_HL6=c(blue2red(5)[-1],"black",rep("grey",ncol(spp[["batt_HL6"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(spp[["upinn"]]))),
  smuti=c(blue2red(5)[c(1,3,4)],rep("grey",ncol(spp[["smuti"]]))),
  cutl=c(blue2red(5)[-2],"black",rep("grey",ncol(spp[["cutl"]]))),
  cera=c(blue2red(5)[-c(2,4)],"black",rep("grey",ncol(spp[["cera"]]))),
  rudi=c(blue2red(5)[-2],"black",rep("grey",ncol(spp[["rudi"]]))),
  ulva=c(blue2red(5),"black",rep("grey",ncol(spp[["ulva"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(spp))
{
circos.clear()
circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)

chordDiagram(x = spp[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(spp[[i]]))))))
mtext(names(spp)[i],line=-5,cex=5)
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
}
dev.off()


##### western North America ####

wNA <- spp[c("giga","gverm","batt","batt_HL1","batt_HL6","upinn","smuti","cutl","cera","rudi","ulva")]
wNA[["giga"]] <- wNA[["giga"]][,c("soCal","NW.America")]
wNA[["gverm"]] <- wNA[["gverm"]][,c("Cali","PNW")]
wNA[["upinn"]] <- as.matrix(data.frame(NAm_south=wNA[["upinn"]][,c("NAm_south")])) 
wNA[["smuti"]] <- wNA[["smuti"]][,c("Cali","PNW")]; 
wNA[["rudi"]] <- as.matrix(data.frame(NAm_north=wNA[["rudi"]][,c("NAm_north")]))
wNA[["ulva"]] <- as.matrix(data.frame(NAm_south=wNA[["ulva"]][,c("NAm_south")]))

cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(wNA[["giga"]]))),
  gverm=c(blue2red(5),"black",rep("grey",ncol(wNA[["gverm"]]))),
  batt=c(blue2red(5)[-1],"black",rep("grey",ncol(wNA[["batt"]]))),
  batt_HL1=c(blue2red(5)[-1],"black",rep("grey",ncol(wNA[["batt_HL1"]]))),
  batt_HL6=c(blue2red(5)[-1],"black",rep("grey",ncol(wNA[["batt_HL6"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(wNA[["upinn"]]))),
  smuti=c(blue2red(5)[c(1,3,4)],rep("grey",ncol(wNA[["smuti"]]))),
  cutl=c(blue2red(5)[-2],"black",rep("grey",ncol(wNA[["cutl"]]))),
  cera=c(blue2red(5)[-c(2,4)],"black",rep("grey",ncol(wNA[["cera"]]))),
  rudi=c(blue2red(5)[-2],"black",rep("grey",ncol(wNA[["rudi"]]))),
  ulva=c(blue2red(5),"black",rep("grey",ncol(wNA[["ulva"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots_wNA.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(wNA))
{
  circos.clear()
  circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)
  
  chordDiagram(x = wNA[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(wNA[[i]]))))))
  mtext(names(wNA)[i],line=-5,cex=5)
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
}
dev.off()

### California ###
cal <- spp[c("giga","gverm","batt","batt_HL1","batt_HL6","upinn","smuti","cutl","ulva")]
cal[["giga"]] <- as.matrix(data.frame(Cali=cal[["giga"]][,c("soCal")]))
cal[["gverm"]] <- as.matrix(data.frame(Cali=cal[["gverm"]][,c("Cali")])) 
cal[["batt"]] <- as.matrix(data.frame(Cali=cal[["batt"]][,c("NAm_south")]))
cal[["batt_HL1"]] <- as.matrix(data.frame(Cali=cal[["batt_HL1"]][,c("Cali")]))
cal[["batt_HL6"]] <- as.matrix(data.frame(Cali=cal[["batt_HL6"]][,c("Cali")]))
cal[["upinn"]] <- as.matrix(data.frame(Cali=cal[["upinn"]][,c("NAm_south")])) 
cal[["smuti"]] <- as.matrix(data.frame(Cali=cal[["smuti"]][,c("Cali")]))
cal[["ulva"]] <- as.matrix(data.frame(Cali=cal[["ulva"]][,c("NAm_south")]))

cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(cal[["giga"]]))),
  gverm=c(blue2red(5),"black",rep("grey",ncol(cal[["gverm"]]))),
  batt=c(blue2red(5)[-1],"black",rep("grey",ncol(cal[["batt"]]))),
  batt_HL1=c(blue2red(5)[-1],"black",rep("grey",ncol(cal[["batt_HL1"]]))),
  batt_HL6=c(blue2red(5)[-1],"black",rep("grey",ncol(cal[["batt_HL6"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(cal[["upinn"]]))),
  smuti=c(blue2red(5)[c(1,3,4)],rep("grey",ncol(cal[["smuti"]]))),
  cutl=c(blue2red(5)[-2],"black",rep("grey",ncol(cal[["cutl"]]))),
  ulva=c(blue2red(5),"black",rep("grey",ncol(cal[["ulva"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots_NAm_south.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(cal))
{
  circos.clear()
  circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)
  
  chordDiagram(x = cal[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(cal[[i]]))))))
  mtext(names(cal)[i],line=-5,cex=5)
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
}
dev.off()

#### Pacific North America ####
pnw <- spp[c("giga","gverm","batt","batt_HL1","batt_HL6","smuti","cera","rudi")]
pnw[["giga"]] <- as.matrix(data.frame(pnwi=pnw[["giga"]][,c("NW.America")]))
pnw[["gverm"]] <- as.matrix(data.frame(pnwi=pnw[["gverm"]][,c("PNW")])) 
pnw[["batt"]] <- as.matrix(data.frame(pnwi=pnw[["batt"]][,c("NAm_north")]))
pnw[["batt_HL1"]] <- as.matrix(data.frame(pnwi=pnw[["batt_HL1"]][,c("PNW")]))
pnw[["batt_HL6"]] <- as.matrix(data.frame(pnwi=pnw[["batt_HL6"]][,c("PNW")]))
pnw[["smuti"]] <- as.matrix(data.frame(pnwi=pnw[["smuti"]][,c("PNW")]))
pnw[["cera"]] <- as.matrix(data.frame(pnwi=pnw[["cera"]][,c("NAm_north")]))
pnw[["rudi"]] <- as.matrix(data.frame(pnwi=pnw[["rudi"]][,c("NAm_north")]))

cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(pnw[["giga"]]))),
  gverm=c(blue2red(5),"black",rep("grey",ncol(pnw[["gverm"]]))),
  batt=c(blue2red(5)[-1],"black",rep("grey",ncol(pnw[["batt"]]))),
  batt_HL1=c(blue2red(5)[-1],"black",rep("grey",ncol(pnw[["batt_HL1"]]))),
  batt_HL6=c(blue2red(5)[-1],"black",rep("grey",ncol(pnw[["batt_HL6"]]))),
  smuti=c(blue2red(5)[c(1,3,4)],rep("grey",ncol(pnw[["smuti"]]))),
  cera=c(blue2red(5)[-c(2,4)],"black",rep("grey",ncol(pnw[["cera"]]))),
  rudi=c(blue2red(5)[-2],"black",rep("grey",ncol(pnw[["rudi"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots_NAm_north.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(pnw))
{
  circos.clear()
  circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)
  
  chordDiagram(x = pnw[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(pnw[[i]]))))))
  mtext(names(pnw)[i],line=-5,cex=5)
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
}
dev.off()


##### Europe ####

Eur <- spp[c("giga","gverm","upinn","smuti","rudi","ulva")]
Eur[["giga"]] <- Eur[["giga"]][,c("EuropeNorth","EuropeSouth")]
Eur[["gverm"]] <- Eur[["gverm"]][,c("EuropeNorth","EuropeSouth")]
Eur[["upinn"]] <- Eur[["gverm"]][,c("EuropeNorth","EuropeSouth")]
Eur[["smuti"]] <- Eur[["smuti"]][,c("EuropeNorth","EuropeSouth")]; Eur[["smuti"]] <- Eur[["smuti"]]+0.09 # to see it.
Eur[["rudi"]] <- as.matrix(data.frame(NAm_north=Eur[["rudi"]][,c("EuropeSouth")]))
Eur[["ulva"]] <- Eur[["ulva"]][,c("EuropeNorth","EuropeSouth")]

cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(Eur[["giga"]]))),
  gverm=c(blue2red(5),"black",rep("grey",ncol(Eur[["gverm"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(Eur[["upinn"]]))),
  smuti=c(blue2red(5)[c(1,3,4)],rep("grey",ncol(Eur[["smuti"]]))),
  rudi=c(blue2red(5)[-2],"black",rep("grey",ncol(Eur[["rudi"]]))),
  ulva=c(blue2red(5),"black",rep("grey",ncol(Eur[["ulva"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots_Eur.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(Eur))
{
  circos.clear()
  circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)
  
  chordDiagram(x = Eur[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(Eur[[i]]))))))
  mtext(names(Eur)[i],line=-5,cex=5)
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
}
dev.off()

##### New Zealand ####

NZ <- spp[c("giga","upinn","ulva")]
NZ[["giga"]] <- as.matrix(data.frame(NZ=NZ[["giga"]][,"New.Zealand"]))
NZ[["upinn"]] <- as.matrix(data.frame(NZ=NZ[["upinn"]][,"New.Zealand"]))
NZ[["ulva"]] <- as.matrix(data.frame(NZ=NZ[["ulva"]][,"NZ"]))

cols.to.use <- list(
  giga=c(blue2red(5),"black",rep("grey",ncol(NZ[["giga"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(NZ[["upinn"]]))),
  upinn=c(blue2red(5),"black",rep("grey",ncol(NZ[["ulva"]]))))

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestCirclizePlots_NZ.pdf",height=30,width=30)
par(mfrow=c(3,3))
for (i in 1:length(NZ))
{
  circos.clear()
  circos.par(gap.after = c(rep(5,5),15,rep(5,5),15),start.degree = 90, gap.degree = 4)
  
  chordDiagram(x = NZ[[i]], directional = 1,grid.col = cols.to.use[[i]], annotationTrack = "grid",transparency = 0.25,annotationTrackHeight = c(0.1, 0.1),self.link=10,direction.type = c("arrows", "diffHeight"),link.arr.type = "big.arrow",link.arr.length =  0.15,diffHeight = -0.001,preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(NZ[[i]]))))))
  mtext(names(NZ)[i],line=-5,cex=5)
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
}
dev.off()

