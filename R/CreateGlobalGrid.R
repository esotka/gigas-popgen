##########################################################
### 1) create global 1º x 1º grid
### borrowed code from Palacios https://rpubs.com/danielequs/199150
##########################################################

rm(list=ls())
library(spatstat) # marks, ppp, quadratcount
library(rgdal) #readOGR
library(tools) #file_path_sans_ex
library(ggplot2) # fortify
library(dplyr) # inner_join
#library(scales) # alpha
#library(utils) 
#
#library(readxl)

path.ne.coast <- ("data/ne_10m_coastline/") # from Natural Earth
fnam.ne.coast <- "ne_10m_coastline.shp"
dat.coast <- readOGR(dsn = path.ne.coast,layer = file_path_sans_ext(fnam.ne.coast))

# Provide the function quick.subset() from Simon Goring's page:
# https://downwithtime.wordpress.com/tag/maps/
quick.subset <- function(x, domain){
  x@data$id <- rownames(x@data)
  x.f = fortify(x, region = "id")
  x.join <- inner_join(x.f, x@data, by = "id")
  x.subset <- subset(x.join, x.join$long > domain[1] & 
                       x.join$long < domain[2] & 
                       x.join$lat > domain[3] & 
                       x.join$lat < domain[4])
  x.subset
}
# domain should be a vector of four values: c(xmin, xmax, ymin, ymax)

domain <- c(-130,180,-50,75) 
# Extract the coastline data for the desired domain using quick.subset():
df.coastline <- quick.subset(dat.coast, domain) 
ppp.coastline <- ppp(x=df.coastline$long,y=df.coastline$lat,xrange=c(-130,180),yrange=c(-50,75))
grid.globe <- quadratcount(ppp.coastline, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4])) ### 1º x 1º latitude blocks; Divides window into quadrats and counts the numbers of points in each quadrat.
lon.to.use <- seq(-129.5,179.5,1) # center of the columns of quadratcount grid
lat.to.use <- sort(seq(-49.5,74.5,1),decreasing=T)
colnames(grid.globe) <- lon.to.use
rownames(grid.globe) <- lat.to.use
grid.globe01 <- ifelse(grid.globe==0,0,1); image(t(grid.globe01)) # show an image of the globe
df.globe <- data.frame(grid.globe) # turn quadrat count into dataframe
df.globe$y <- as.numeric(as.character(df.globe$y))
df.globe$x <- as.numeric(as.character(df.globe$x))
colnames(df.globe)[1:2] <- c("lat","lon")
df.globe$coastline <- ifelse(df.globe$Freq==0,0,1)

### here's the set of coastlines to model. Chosen based on distribution of gigas
grids.to.use <- list(
  asia_domain = c(105,145,30,45),
  wNA_domain = c(-130,-111,23,60),
  domain.eur1 = c(-12.14,-5.27,31.87,44.65), # eur: only Atlantic Ocean pops. Ignoring Mediterranean for now.
  domain.eur2 = c(-12.14,5.46,44.65,59.56),
  domain.eur3 = c(5.46,14.90,52.91,59.56),
  domain.eur4 = c(-6.96,0.76,43.06,45.06))

out <- list()
for (j in 1:length(grids.to.use))
{
  tmp <- grids.to.use[[j]]
  coastline <- quick.subset(dat.coast, tmp) 
  ppp.tmp <- ppp(x=coastline$long,y=coastline$lat,xrange=c(-130,180),yrange=c(-50,75))
  counts.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  out[[names(grids.to.use)[j]]] <- data.frame(counts.tmp)$Freq
}

out2 <- ifelse(rowSums(matrix(unlist(out),nrow=38750))>0,1,0)
df.globe$grid.for.model <- out2

# name each quadrat
df.globe$reg <- c()
df.globe$reg[df.globe$lon <(-100)] <- "2_wNA" 
df.globe$reg[df.globe$lon < 50 & df.globe$lon > (-50)] <- "3_Eur" 
df.globe$reg[df.globe$lon >100] <- "1_Asia" 
df.globe$gridID <- 1:dim(df.globe)[1]


# print map of grid for model
df.coastline <- df.globe[df.globe$coastline==1,]
df.onlyModelgrid <- df.coastline[df.coastline$grid.for.model==1,] # red points are model
plot(x=df.coastline$lon,y=df.coastline$lat,cex=.01,pch=20,xlab="",ylab="",main="red points = grid for model")
points(x=df.onlyModelgrid$lon,y=df.onlyModelgrid$lat,col="red",pch=20,cex=0.1)

plot(x=df.onlyModelgrid$lon,y=df.onlyModelgrid$lat,cex=.1,pch=20,xlab="",ylab="",col="red",main="red points = grid for model")

# make the modelGrid into a square matrix (22x22)
df.onlyModelgrid[448:462,] <- rep(NA,ncol(df.onlyModelgrid))
df.onlyModelgrid$reg <- factor(df.onlyModelgrid$reg)
df.onlyModelgrid <- df.onlyModelgrid[order(df.onlyModelgrid$reg),]
df.onlyModelgrid$plotted <- expand.grid(1:22,1:21)

plot(df.onlyModelgrid$plotted$Var1,df.onlyModelgrid$plotted$Var2,pch=21,col=alpha(c("black","red","blue")[df.onlyModelgrid$reg],.4),yaxt="n",xaxt="n",ylab="",xlab="",main="Asia (black); wNA (red); Europe (blue) " )

## write grid for empirical data
df.globe$plotted.Var1 <- df.onlyModelgrid$plotted$Var1[match(df.globe$gridID,df.onlyModelgrid$gridID)]
df.globe$plotted.Var2 <- df.onlyModelgrid$plotted$Var2[match(df.globe$gridID,df.onlyModelgrid$gridID)]
df.globe <- df.globe[,-3]
save(df.globe,file="output/df.globe.Rda") 

