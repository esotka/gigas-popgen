### make a map with regions

rm(list=ls())
#library(readxl,warn.conflicts = F,quietly = T)
#library(spatstat,warn.conflicts = F,quietly = T) # marks, ppp, quadratcount
library(scales,warn.conflicts = F,quietly = T) # alpha
#library(circlize)
library(colorRamps)
#library(reshape)
library(maps)
library(mapdata)
#library(ranger)


domain <- c(-130,180,-50,75) 
load("../df.globe.Rda")
meta_source <- read.csv("../df.globe_source.csv")

map1 <- function() {
  map("worldHires",xlim=c(125,155),ylim=c(28,47),col="gainsboro",fill=TRUE)
  tmp <- meta_source[!meta_source$sourceID=="nonSource",]
  tmp$sourceID <- factor(tmp$sourceID)
  tmp$sourceID <- factor(tmp$sourceID,levels=levels(tmp$sourceID)[c(1,2,5,4,3)])
  rect(tmp$lon.sq1,tmp$lat.sq2,tmp$lon.sq2,tmp$lat.sq1,col=alpha(blue2red(5),.5)[tmp$sourceID])
  segments(x0=seq(130,155,10),y0=27,x1=seq(130,155,10),y1=28.3)
  text(x=c(130,140,150),y=28.7,paste(c(130,140,150),"º"),cex=.8)
  segments(124,seq(30,45,5),125.5,seq(30,45,5))
  text(x=125,y=seq(30,45,5),paste(seq(30,45,5),"º"),cex=.8,pos=4)
  text(x=c(144.7672,143.0341,140.8678,134.6165,132), y=c(41.37038, 38.15945, 34.36471,32,29.5),c("Hokkaido","Miyagi","Tokyo","Seto Inland Sea","Kagoshima"),col=blue2red(5),pos=4,cex=1.2)
  box()
#  points(meta$longitude,meta$latitude,pch=20,cex=2)
}
png("JapanRegions_Map.png",height=5,width=5,units="in",res=400); map1(); dev.off()
