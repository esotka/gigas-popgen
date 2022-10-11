load("FINALFIGS/0_globalGrid/df.globe.Rda")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
reg.to.use <- matrix(c(
  "Japan","Native",
  "soCal","wNA",
  "France","Eur",
  "Spain","Eur",
  "Korea","Native",
  "NW America","wNA",
  "Ireland","Eur",
  "Sweden","Eur",
  "Norway","Eur",
  "Denmark","Eur",
  "Europe","Eur"),ncol=2,byrow = T)
tmp <- meta[meta$Region%in%c(reg.to.use[,1]),]
tmp$Region2 <- reg.to.use[match(tmp$Region,reg.to.use[,1]),2]
## find which quadrat each pop
gridIDs.tmp <- c()
for (i in 1:length(tmp$pop))
{
  ppp.tmp <- ppp(x=tmp$Longitude[tmp$pop==tmp$pop[i]],y=tmp$Latitude[tmp$pop==tmp$pop[i]],xrange=domain[1:2],yrange=domain[3:4])
  grid.tmp <- quadratcount(ppp.tmp, nx = abs(domain[1]-domain[2]), ny = abs(domain[3]-domain[4]))
  df.tmp <- data.frame(grid.tmp,quadIDs=1:dim(df.globe)[1])
  gridIDs.tmp <- c(gridIDs.tmp,df.tmp$quadIDs[df.tmp$Freq==1])
}
gridIDs.meta <- data.frame(gridIDs=gridIDs.tmp,pop=tmp$pop)
gridIDs.meta$vars <- df.globe[match(gridIDs.meta$gridIDs,df.globe$gridID),]

plot(df.globe$plotted.Var1,df.globe$plotted.Var2,pch=21,col=alpha(c("black","red","blue")[factor(df.globe$reg)],.4),yaxt="n",xaxt="n",ylab="",xlab="",main="Asia (black); wNA (red); Europe (blue) " )

points(gridIDs.meta$vars$plotted.Var1,gridIDs.meta$vars$plotted.Var2,pch=20,cex=2,col=alpha(c("black","red","blue")[as.factor(gridIDs.meta$vars$reg)],.4))
write.csv(gridIDs.meta,"FINALFIGS/5_circlize/cgiga/cgiga_meta.csv",quote=F,row.names = F)
