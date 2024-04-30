### make pretty SST plot
### SST + SAT from BioOracle marine dataset.
rm(list=ls())
library(raster)
library(gplots)
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
cols.to.use <- c(blue2red(5))
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima")
meta$cols.to.use = cols.to.use[match(meta$Region2,names(cols.to.use))]
meta$cols.to.use[is.na(meta$cols.to.use)] = "grey"
meta$cols.to.use[meta$Region=="Korea"] = "black"
meta$popTextCol <- c("black","white")[factor(meta$Region2=="nonSource")]
meta$popTextCol[meta$Region2=="Akkeshi"] <- "white"

Grids <- raster::stack("FINALFIGS/8_mapSamples/sstmean.asc")
nat.range <- c(110,160,25,50)
Grids.nat <- crop(Grids,nat.range)
rdf <- as.data.frame(Grids.nat, xy=TRUE)
g1 <- ggplot(data = rdf) +
            geom_raster(mapping=aes(x=x, y=y, fill=sstmean))+ scale_fill_gradientn(colours= alpha(bluered(100),.5)) +
            xlim(123,160) + ylim(25,47) + xlab("") + ylab("") +
            annotate("segment", x = meta$Longitude, xend = meta$lon.pie, y = meta$Latitude, yend = meta$lat.pie) +
            geom_point(data = meta, aes(x=lon.pie,y=lat.pie), col=meta$cols.to.use,cex=6.5) +
            geom_point(data = meta, aes(x=lon.pie,y=lat.pie), col="black",cex=6.5,pch=21) + 
            annotate("text",x=meta$lon.pie,y=meta$lat.pie,label=meta$pop,colour=meta$popTextCol,cex=2) +
            theme(panel.background = element_blank())

png('FINALFIGS/8_mapSamples/japan_SST.png',width=7,height=4,units="in",res=500)
print(g1)
dev.off()

