### points of PC1 output of SNPs
rm(list=ls())
library(maps);library(mapdata);library(colorRamps)#library(fields)#;library(reshape);;library(plotrix);library(scales)
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")

cols.to.use <- c(blue2red(5))
names(cols.to.use) <- c("Akkeshi","Miyagi","Tokyo","Seto","Kagoshima")
meta$cols.to.use = cols.to.use[match(meta$Region2,names(cols.to.use))]
meta$cols.to.use[is.na(meta$cols.to.use)] = "grey"
meta$cols.to.use[meta$Region=="Korea"] = "black"
meta$popTextCol <- c("black","white")[factor(meta$Region2=="nonSource")]
meta$popTextCol[meta$Region2=="Akkeshi"] <- "white"

### Japan
#quartz(width=5,height=4)
png('FINALFIGS/8_mapSamples/japan.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map("worldHires","China",bg="lightsteelblue1",fill=T,col="white",xlim=c(123,150),ylim=c(28,47),xaxt="n",yaxt="n",ylab="",xlab="")
map("worldHires","South Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","North Korea",add=TRUE,col="white",fill=TRUE)
map("worldHires","Japan",add=TRUE,col="white",fill=TRUE)
map("worldHires","USSR",col="white",fill=TRUE,add=TRUE)
segments(seq(125,145,5),25,seq(125,145,5),28.5)
text(x=c(125,145),y=28,paste(c(125,145),"º",sep=""),pos=3,cex=.5)
segments(120,seq(30,45,5),124,seq(30,45,5))
text(x=123.5,c(30,45)-0.15,paste(c(30,45),"º",sep=""),pos=4,cex=.5)
box()
segments(meta$Longitude,meta$Latitude,meta$lon.pie,meta$lat.pie,col="black",lwd=1)
points(x=meta$lon.pie,y=meta$lat.pie,col=meta$cols.to.use,pch=20,cex=4)
points(x=meta$lon.pie,y=meta$lat.pie,col="black",pch=21,cex=2.8)
text(x=meta$lon.pie,y=meta$lat.pie,meta$pop,col=meta$popTextCol,cex=.5)
dev.off()


### Europe range
png('FINALFIGS/8_mapSamples/europe.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map('worldHires',xlim=c(-17,24),ylim=c(34,67),bg="lightsteelblue1",fill=T,col="white")
box()
segments(meta$Longitude,meta$Latitude,meta$lon.pie,meta$lat.pie,col="black",lwd=1)
points(x=meta$lon.pie,y=meta$lat.pie,col=meta$cols.to.use,pch=20,cex=4)
points(x=meta$lon.pie,y=meta$lat.pie,col="black",pch=21,cex=2.8)
text(x=meta$lon.pie,y=meta$lat.pie,meta$pop,col="black",cex=.5)
dev.off()

## pac nw
#quartz(width=5,height=4)
png('FINALFIGS/8_mapSamples/PNW.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map("worldHires",xlim=c(-125,-120),ylim=c(46,50),bg="lightsteelblue1",fill=TRUE,col="white")
points(meta$Longitude, meta$Latitude, pch=19, col="grey", cex=2.9)  
points(x=meta$Longitude,y=meta$Latitude,col="black",pch=21,cex=2.8)
text(meta$Longitude,meta$Latitude,meta$pop,cex=.5,col="black")
box()
segments(seq(-124,-121,1),45,seq(-124,-121,1),46.2)
text(x=seq(-123,-121,1),y=46.1,paste(seq(-123,-121,1),"º",sep=""),pos=3,cex=.5)
segments(-120.2,seq(44,52,1),-117,seq(44,52,1))
text(x=-120.7,y=c(47,48),paste(c(47,48),"º",sep=""),pos=4,cex=.5)
text(x=-121.5,y=49.5,"Canada")
text(x=-121,y=48.5,"USA (WA)",cex=.8)
dev.off()

# so Cal
#quartz(width=5,height=4)
png('FINALFIGS/8_mapSamples/soCal.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map("worldHires",xlim=c(-120,-115),ylim=c(31,35),bg="lightsteelblue1",fill=TRUE,col="white")
tmp = meta[meta$pop%in%c("ALB","NEW","PES"),]
points(tmp$Longitude, tmp$Latitude, pch=19, col="grey", cex=2.9)
points(x=tmp$Longitude,y=tmp$Latitude,col="black",pch=21,cex=2.8)
text(tmp$Longitude,tmp$Latitude,tmp$pop,cex=.5,col="black")
tmp = meta[meta$pop%in%c("GRC"),]
segments(tmp$Longitude-.5,tmp$Latitude,tmp$Longitude,tmp$Latitude,col="black")
points(tmp$Longitude-.5, tmp$Latitude, pch=19, col="grey", cex=3)
points(x=tmp$Longitude-.5,y=tmp$Latitude,col="black",pch=21,cex=2.8)
text(tmp$Longitude-.5,tmp$Latitude,tmp$pop,cex=.5,col="black")
tmp = meta[meta$pop%in%c("TJE"),]
segments(tmp$Longitude-.5,tmp$Latitude-.4,tmp$Longitude,tmp$Latitude,col="black")
points(tmp$Longitude-.5, tmp$Latitude-.4, pch=19, col="grey", cex=3)
points(x=tmp$Longitude-.5,y=tmp$Latitude-.4,col="black",pch=21,cex=2.8)
text(tmp$Longitude-.5,tmp$Latitude-.4,tmp$pop,cex=.5,col="black")
box()
segments(seq(-120,-115,1),30,seq(-120,-115,1),31.2)
text(x=seq(-120,-115,1),y=31.1,paste(seq(-120,-115,1),"º",sep=""),pos=3,cex=.5)
segments(-121,seq(30,35,1),-119.7,seq(30,35,1))
text(x=-119.7,y=seq(30,35,1),paste(seq(30,35,1),"º",sep=""),pos=4,cex=.5)
text(x=-116,y=33.5,"California")
dev.off()

# NZ
#quartz(width=5,height=4)
png('FINALFIGS/8_mapSamples/NZ.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map("worldHires",xlim=c(170,180),ylim=c(-42,-35),bg="lightsteelblue1",fill=TRUE,col="white")
points(meta$Longitude, meta$Latitude, pch=19, col="grey", cex=2.9)
points(x=meta$Longitude,y=meta$Latitude,col="black",pch=21,cex=2.8)
text(meta$Longitude,meta$Latitude,meta$pop,cex=.5,col="black")
box()
segments(seq(170,180,2),-44,seq(170,180,2),-41.5)
text(x=seq(170,180,2),y=-41.5,paste(seq(170,180,2),"º",sep=""),pos=3,cex=.5)
segments(170.5,seq(-42,-35,2),169,seq(-42,-35,2))
text(x=170.5,y=seq(-42,-35,2),paste(seq(-42,-35,2),"º",sep=""),pos=4,cex=.5)
dev.off()

# Argentina
#quartz(width=5,height=4)
png('FINALFIGS/8_mapSamples/Argentina.png',width=5,height=4,units="in",res=500)
par(mar=c(0,0,0,0))
map("worldHires",xlim=c(-70,-55),ylim=c(-50,-35),bg="lightsteelblue1",fill=TRUE,col="white")
points(meta$Longitude, meta$Latitude, pch=19, col="grey", cex=3)
points(x=meta$Longitude,y=meta$Latitude,col="black",pch=21,cex=2.8)
text(meta$Longitude,meta$Latitude,meta$pop,cex=.5,col="black")
box()
segments(seq(-70,-55,5),-51,seq(-70,-55,5),-49)
text(x=seq(-70,-55,5),y=-49,paste(seq(-70,-55,5),"º",sep=""),pos=3,cex=.5)
segments(-72,seq(-50,-35,5),-69,seq(-50,-35,5))
text(x=-69.5,y=seq(-50,-35,5),paste(seq(-50,-35,5),"º",sep=""),pos=4,cex=.5)
dev.off()


### globe
#quartz(width=10,height=8)
png('FINALFIGS/8_mapSamples/globe.png',width=10,height=8,units="in",res=500)
map('worldHires')
tmp <- meta[meta$Region2%in%c("Argentina","Chile","New Zealand","PNW","soCalifornia"),]
lat <- tapply(tmp$Latitude,tmp$Region2,mean)
lon <- tapply(tmp$Longitude,tmp$Region2,mean)
#cols.to.use <- tmp$pc1.cols[match(names(lat),tmp$GeneticRegions)]
points(x=lon,y=lat,col="grey",pch=20,cex=4)
points(x=lon,y=lat,col="black",pch=21,cex=3)
dev.off()

