### pie diagram with all the species
# make a PDF with all chordDiagrams
# make sure you run all of the species-specific RF files first (e.g., rudiRF&Circlize.R)
#library(circlize)
library(colorRamps)
rm(list=ls())

filenames <- 
  c("cgigasByReg.csv",
    "upinnByReg.csv",
    "dideByReg.csv",
    "gvermSNPByReg.csv",
    "hamiByReg.csv",
    "battrByReg.csv",
    "battr_HL1ByReg.csv",
    "battr_HL6ByReg.csv",
    "mcylByReg.csv",
    "takaByReg.csv",
    "sangByReg.csv")

sppnames <- sub("ByReg.csv","",x=filenames)

spp <- list()
for (i in 1:length(filenames))
{
  spp[[i]] <- read.csv(paste("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/",filenames[i],sep=""))  
  rownames(spp[[i]]) <- spp[[i]][,1]
  spp[[i]] <- as.matrix(data.frame(spp[[i]][-1]))
#  if(sppnames[i]%in%c("mcyl","taka","sang") {spp[[i]] <- as.matrix(data.frame(spp[[i]][-1]))}
#  else{spp[[i]] <- as.matrix(spp[[i]][,-1])}
}
names(spp) <- sppnames

pdf("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/RandomForestPieCharts.pdf",width=6,height=12)
par(mfrow=c(11,5),mar=c(0,0,0,1))
#cgigas
tmp <- spp[["cgigas"]]
col.to.use <- c(blue2red(5),"black")
pie(tmp[,c("NW.America")],col = col.to.use,labels = "",main="nAM_north") #nAM_north
pie(tmp[,c("soCal")],col = col.to.use,labels = "",main="nAM_south") #nAM_south
pie(tmp[,c("EuropeNorth")],col = col.to.use,labels = "",main="Europe_north") #EuropeNorth
pie(tmp[,c("EuropeSouth")],col = col.to.use,labels = "",main="Europe_south") #EuropeSouth
pie(tmp[,c("New.Zealand")],col = col.to.use,labels = "",main="NZ") #NZ

#upinn
tmp <- spp[["upinn"]]
col.to.use <- c(blue2red(5),"black")
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #nAM_north
pie(tmp[,c("USA")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #EuropeNorth
pie(tmp[,c("France")],col = col.to.use,labels = "") #EuropeSouth
pie(tmp[,c("New.Zealand")],col = col.to.use,labels = "") #NZ

#dide
tmp <- spp[["dide"]]
col.to.use <- c(blue2red(5)[1:3])#,"black")
pie(tmp[,c("NAm_north")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("NAm_south")],col = col.to.use,labels = "") #nAM_south
pie(tmp[,c("EuropeNorth")],col = col.to.use,labels = "") #EuropeNorth
pie(tmp[,c("EuropeSouth")],col = col.to.use,labels = "") #EuropeSouth
pie(tmp[,c("NewZealand")],col = col.to.use,labels = "") #NZ

#gverm
tmp <- spp[["gvermSNP"]]
col.to.use <- c(blue2red(5)[-4],"black")
pie(tmp[,c("WNA")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("soCal")],col = col.to.use,labels = "") #nAM_south
pie(tmp[,c("EuropeNorth")],col = col.to.use,labels = "") #EuropeNorth
pie(tmp[,c("EuropeSouth")],col = col.to.use,labels = "") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#hami
tmp <- spp[["hami"]]
col.to.use <- c(blue2red(5),"black")
pie(tmp[,c("Washington")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("California")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #EuropeNorth
pie(rowSums(tmp[,c("France","Italy","Spain")]),col = col.to.use,labels = "") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#battr
tmp <- spp[["battr"]]
col.to.use <- c(blue2red(5)[-1],"black")
pie(tmp[,c("NAm_north")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("NAm_south")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeNorth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#battr_HL1
tmp <- spp[["battr_HL1"]]
col.to.use <- c(blue2red(5)[-1],"black")
pie(tmp[,c("NAm_north")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("NAm_south")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeNorth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#battr_HL6
tmp <- spp[["battr_HL6"]]
col.to.use <- c(blue2red(5)[-1],"black")
pie(tmp[,c("NAm_north")],col = col.to.use,labels = "") #nAM_north
pie(tmp[,c("NAm_south")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeNorth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#mcyl
tmp <- spp[["mcyl"]]
col.to.use <- c(blue2red(5),"black")
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #nAM_north
pie(tmp[,c("USA")],col = col.to.use,labels = "") #nAM_south
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeNorth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#taka
tmp <- spp[["taka"]]
col.to.use <- c(blue2red(5),"black")
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #nAM_north
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #nAM_south
pie(tmp[,c("NorthEurope")],col = col.to.use,labels = "") #EuropeNorth
pie(tmp[,c("SouthEurope")],col = col.to.use,labels = "") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ

#sang
tmp <- spp[["sang"]]
col.to.use <- c(blue2red(5),"black")
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #nAM_north
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,cex=3,col="grey","NA") #nAM_south
pie(tmp[,c("Europe")],col = col.to.use,labels = "") #EuropeNorth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #EuropeSouth
pie(1,col = "white",border="white",labels = ""); mtext(line=-6,side=3,col="grey",cex=3,"NA") #NZ


dev.off()

