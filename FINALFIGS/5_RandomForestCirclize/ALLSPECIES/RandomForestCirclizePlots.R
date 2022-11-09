# make a PDF with all chordDiagrams
# make sure you run all of the species-specific RF files first (e.g., rudiRF&Circlize.R)
library(circlize)
library(colorRamps)
rm(list=ls())

filenames <- c(
  "cgigasByReg.csv",
  "gvermSNPByReg.csv",
  "upinnByReg.csv",
  "mcylByReg.csv",
  "battr_HL6ByReg.csv",
  "battr_HL1ByReg.csv",
  "battrByReg.csv")
sppnames <- c(
  "Crassostrea_gigas",
  "Gracilaria_vermiculophylla",
  "Undaria_pinnitifida",
  "Mutimo_cylindricus",
  "Trematode_HL1",
  "Trematode_HL6",
  "Batillaria_attramentria"
)
spp <- list()
for (i in 1:length(filenames))
{
spp[[i]] <- read.csv(paste("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/",filenames[i],sep=""))  
rownames(spp[[i]]) <- spp[[i]][,1]; 
if(sppnames[i]=="Mutimo_cylindricus"){spp[[i]] <- as.matrix(data.frame(spp[[i]][-1]))}
  else{spp[[i]] <- as.matrix(spp[[i]][,-1])}
}
names(spp) <- sppnames

cols.to.use <- list(
  c(blue2red(5),"black",rep("grey",ncol(spp[[1]]))),# Crassostrea_gigas
  c(blue2red(5)[-4],"black",rep("grey",ncol(spp[[2]]))),# Gracilaria_vermiculophylla
  c(blue2red(5),"black",rep("grey",ncol(spp[[3]]))),# Undaria_pinnitifida
  c(blue2red(5),"black",rep("grey",1)),# Mutimo_cylindricus
  c(blue2red(5)[-1],"black",rep("grey",ncol(spp[[5]]))),# Trematode_HL1
  c(blue2red(5)[-1],"black",rep("grey",ncol(spp[[6]]))),# Trematode_HL6
  c(blue2red(5)[-1],"black",rep("grey",ncol(spp[[7]]))))
# Batillaria_attramentria

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

wNA <- list()
wNA[["Crassostrea_gigas"]] = spp[["Crassostrea_gigas"]][,c("soCal","NW.America")]
colnames(wNA[["Crassostrea_gigas"]]) = c("NAm_south","NAm_north")
wNA[["Gracilaria_vermiculophylla"]] = spp[["Gracilaria_vermiculophylla"]][,c("soCal","WNA")]
colnames(wNA[["Gracilaria_vermiculophylla"]]) = c("NAm_south","NAm_north")
wNA[["Undaria_pinnitifida"]] = as.matrix(data.frame(NAm_south=spp[["Undaria_pinnitifida"]][,c("USA")])) 
colnames(wNA[["Undaria_pinnitifida"]]) = "NAm_south"
wNA[["Mutimo_cylindricus"]] <- as.matrix(data.frame(spp[["Mutimo_cylindricus"]][,c("USA")])) 
colnames(wNA[["Mutimo_cylindricus"]]) = "NAm_south"
wNA[["Trematode_HL1"]] = spp[["Trematode_HL1"]]
wNA[["Trematode_HL6"]] = spp[["Trematode_HL6"]]
wNA[["Batillaria_attramentria"]] = spp[["Batillaria_attramentria"]]

cols.to.use <- list(
  c(blue2red(5),"black",rep("grey",2)),# Crassostrea_gigas
  c(blue2red(5)[-4],"black",rep("grey",2)),# Gracilaria_vermiculophylla
  c(blue2red(5),"black",rep("grey",1)),# Undaria_pinnitifida
  c(blue2red(5),"black",rep("grey",1)),# Mutimo_cylindricus
  c(blue2red(5)[-1],"black",rep("grey",2)),# Trematode_HL1
  c(blue2red(5)[-1],"black",rep("grey",2)),# Trematode_HL6
  c(blue2red(5)[-1],"black",rep("grey",2)))
# Batillaria_attramentria


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

##### Europe ####
Eur <- list()
Eur[["Crassostrea_gigas"]] = spp[["Crassostrea_gigas"]][,c("EuropeNorth","EuropeSouth")]
Eur[["Gracilaria_vermiculophylla"]] = spp[["Gracilaria_vermiculophylla"]][,c("EuropeNorth","EuropeSouth")]
Eur[["Undaria_pinnitifida"]] = as.matrix(data.frame(EuropeSouth=spp[["Undaria_pinnitifida"]][,c("France")]))

cols.to.use <- list(
  c(blue2red(5),"black",rep("grey",2)),# Crassostrea_gigas
  c(blue2red(5)[-4],"black",rep("grey",2)),# Gracilaria_vermiculophylla
  c(blue2red(5),"black",rep("grey",1)))# Undaria_pinnitifida
  
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
NZ <- list()
NZ[["Crassostrea_gigas"]] = as.matrix(data.frame(NZ=spp[["Crassostrea_gigas"]][,c("New.Zealand")]))
NZ[["Undaria_pinnitifida"]] = as.matrix(data.frame(NZ=spp[["Undaria_pinnitifida"]][,c("New.Zealand")]))


cols.to.use <- list(
  c(blue2red(5),"black",rep("grey",1)),
  c(blue2red(5),"black",rep("grey",1)))

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

