## get environmental variables

### PCA environment vs PCA genetics
rm(list=ls())
library(sdmpredictors) # pull in data from BioOracle
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
my.sites <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")

my.sites.ssdata <- data.frame(my.sites,extract(ssdata,my.sites[,c("Longitude","Latitude")]))
my.sites.all <- data.frame(my.sites.ssdata,extract(climdata,my.sites[,c("Longitude","Latitude")]))
#my.sites.all <- my.sites.all[complete.cases(my.sites.all),]
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio1"] <- "MeanAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio5"] <- "MaxAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio6"] <- "MinAirTemp"
colnames(my.sites.all)[colnames(my.sites.all)=="WC_bio7"] <- "RangeAirTemp"
