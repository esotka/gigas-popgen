### PCA environment vs PCA genetics

# environmental PCA on populations we sampled


rm(list=ls())
library(sdmpredictors) # pull in data from BioOracle
library(scales) ## transparent points
library(RColorBrewer) ### colors
library(colorRamps) ### colors
library(raster) # extract()
library(GGally) # correlations
library(lattice)
# BIO-ORACLE: 5km resolution
layers.bio2 <- list_layers( datasets="Bio-ORACLE",) 
ssdata <- load_layers(layercodes = c(
  "BO_sstmax", 
  "BO_sstmin", 
  "BO_sstmean",
  "BO_sstrange",
  "BO_salinity",
  "BO_chlomax",
  "BO_chlomin",
  "BO_chlomean",
  "BO_chlorange"
) , equalarea=FALSE, rasterstack=TRUE,datadir = "FINALFIGS/") 
# WorldClim
layers.clim <- list_layers( datasets="WorldClim")
climdata <- load_layers(layercodes=c(
  "WC_bio1", # Annual mean temperature
  "WC_bio5", # Maximum temperature
  "WC_bio6", # Minimum temperature
  "WC_bio7") # Annual temperature range
  #  "WC_bio12", # Annual precipitation
  #  "WC_bio15" # Preciptiation seasonality
  , equalarea=FALSE, rasterstack=TRUE,datadir = "FINALFIGS/")
my.sites <- read.csv("FINALFIGS/PC.out.csv")
my.sites$natnon <- ifelse(my.sites$GeneticRegions%in%c("Akkeshi","Kagoshima","Korea","Miyagi","Seto","Tokyo_wJap"),"Native","Introduced")
#my.sites$lon <- meta$Longitude[match(my.sites$Group.1,meta$pop)]
#my.sites$lat <- meta$Latitude[match(my.sites$Group.1,meta$pop)]
my.sites$Latitude[my.sites$Group.1=="ESN"] <- 46.24123; my.sites$Longitude[my.sites$Group.1=="ESN"] <- -1.219623 # in an estuary; this is the nearest ocean spot

my.sites$Longitude[my.sites$Group.1=="ALB"] <- -118.135526; my.sites$Latitude[my.sites$Group.1=="ALB"] <- 33.746239 # in an estuary; this is the nearest ocean spot
#my.sites$geneticRegion <- meta$GeneticRegions[match(my.sites$Group.1,meta$pop)]

my.sites.ssdata <- data.frame(my.sites,extract(ssdata,my.sites[,c("Longitude","Latitude")]))
my.sites.all <- data.frame(my.sites.ssdata,extract(climdata,my.sites[,c("Longitude","Latitude")]))
#my.sites.all <- my.sites.all[complete.cases(my.sites.all),]
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio1"] <- "MeanAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio5"] <- "MaxAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio6"] <- "MinAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio7"] <- "RangeAirTemp"



## make PCA
pca1 <- prcomp(my.sites.all[,-(1:11)])#,scannf = FALSE, nf = 2)
plot(pca1$x[,1],pca1$x[,2],pch=c(21,19)[factor(my.sites.all$natnon)],col=c("red","black")[factor(my.sites.all$natnon)],xaxt="none",yaxt="none",ylab="",xlab="",cex=3)
text(pca1$x[,1],pca1$x[,2],col=c("red","white")[factor(my.sites.all$natnon)],my.sites.all$Group.1,cex=0.5)
mtext(side=1,line=1,"PC1"); mtext(side=2,line=1,"PC2")
segments(x0=rep(0,13),y0=rep(0,13),x1=pca1$rotation[,1]*20,y1=pca1$rotation[,2]*20)
text(x=pca1$rotation[,1]*20,y=pca1$rotation[,2]*20,rownames(pca1$rotation),col="darkgrey")
legend(x=-32,y=-15,legend=c("Introduced","Native"),pch=c(21,19),col=c("red","black"))




## PC1 loading variable ranks: 1) ChlA max, 2) ChlA range, 3) ChlA mean


print(sort(abs(pca1$rotation[,1]),decreasing = T))


## PC2 loading variable ranks: 1) SSTrange 2) Min Air Temp, 3) Range air temp 4) SST min


print(sort(abs(pca1$rotation[,2]),decreasing = T))


## Correlations among variables with PCs
### Native and Introduced combined
### PC1 best correlates with 
#1) SSTmean (.773)  
#2) Mean Air Temp (.693)  
#3) SSTmin (0.681)  


pdf("FINALFIGS/GenEnvCorrelations.pdf",height=10,width=10)
print(ggscatmat(data=my.sites.all,columns=c(2,12:24),alpha=.1)) # PC1 and biotic  variables
dev.off()
print(sort(cor(my.sites.all[,c(2,12:24)])[1,]))
heatmap(abs(cor(my.sites.all[,c(2,12:24)])),main="native and introduced")



### Native only 
### PC1 best correlates with 
#1) Max air temp (0.787)
#2) SST mean (0.777)
#3) Min Air Temp (0.751)


print(sort(cor(my.sites.all[my.sites.all$natnon=="Native",c(2,12:24)])[1,]))
heatmap(abs(cor(my.sites.all[my.sites.all$natnon=="Native",c(2,12:24)])),main="native only")


### Introduced only 
### PC1 best correlates with 
#1) SST mean (0.744)
#2) SST min (0.714)
#3) Min air temp (0.707)


print(sort(cor(my.sites.all[my.sites.all$natnon=="Introduced",c(2,12:24)])[1,]))
heatmap(abs(cor(my.sites.all[my.sites.all$natnon=="Introduced",c(2,12:24)])),main="Introduced only")



## plot of SSTmean (seen in all three comparisons) that best correlate with PC1


mysettings <- trellis.par.get()
mysettings$strip.background$col <- "lightgrey"
mysettings$box.dot$cex <- 0
trellis.par.set(mysettings)
xyplot(PC1~BO_sstmean | natnon,data=my.sites.all,type=c("p","r"),main="r=0.773 (all); 0.777 (native); 0.744 (introduced)")
print(cor.test(~PC1+BO_sstmean,data=my.sites.all))
print(cor.test(~PC1+BO_sstmean,data=my.sites.all[my.sites.all$natnon=="Native",]))
print(cor.test(~PC1+BO_sstmean,data=my.sites.all[my.sites.all$natnon=="Introduced",]))
write.csv(file="FINALFIGS/PC.out+EnvData.csv",my.sites.all,quote=F,row.names=F)


