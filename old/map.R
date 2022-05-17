### map of sites - globe

library(maps)
library(maptools)
library(scales) ## transparency with alpha
rm(list=ls())
data(wrld_simpl)
meta <- read.csv('data/MetaData_Cgigas.csv')
pdf('output/maps.pdf')
plot(wrld_simpl, col = "grey", axes = T);points(meta$Longitude, meta$Latitude, pch=19, col=alpha("red",.5), cex=1)  
dev.off()