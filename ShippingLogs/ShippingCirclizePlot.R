# Shipping circlize plot
rm(list=ls())
library(circlize)
library(colorRamps)

load("ShippingLogs/shipping_routes.rda")
print(tbl)
mat <- as.matrix(tbl)
mat <- mat[c(1,2,4,3,5),]
mat <- mat+0.1
rownames(mat) <- c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima")
colnames(mat) <- c("soEurope","noEurope","NZ","nAM_north","nAM_south")

cols.to.use <- c(blue2red(5),rep("grey",ncol(mat)))

pdf("ShippingLogs/ShippingCirclizePlot.pdf",width=10,height=10)
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
mtext("shipping",line=-5,cex=2)
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
