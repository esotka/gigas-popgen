##### northwestern North America ####

wNA <- spp[c(1,)]
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
