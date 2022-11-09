### population sizes
library(grid)
library(gridExtra)
library(gtable)
library(reshape)

rm(list=ls())
filenames <- c(
  "cgigas_sampleSize.csv",
  "gvermSNP_sampleSize.csv",
  "upinn_sampleSize.csv",
  "mcyl_sampleSize.csv",
  "battr_HL1_sampleSize.csv",
  "battr_HL6_sampleSize.csv",
  "battr_sampleSize.csv")

spp <- c(
  "Crassostrea_gigas",
  "Gracilaria_vermiculophylla",
  "Undaria_pinnitifida",
  "Mutimo_cylindricus",
  "Trematode_HL1",
  "Trematode_HL6",
  "Batillaria_attramentria"
)
out <- c()
for (i in 1:length(filenames))
{
  out <- rbind(out,data.frame(spp=spp[i],
                              read.csv(paste("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/",filenames[i],sep=""))))
}
colnames(out) <- c("spp","pop","n","Region")
#out$Region[out$Region%in%c("Cali","soCal")] <- "NAm_south"
#out$Region[out$Region%in%c("PNW","NW America")] <- "NAm_north"
NatNon <- ifelse(out$Region%in%c("hok","hon","kag","sea","tok","nonSource"),"Native","Introduced")
out$Region <- paste(NatNon,out$Region,sep="_")


a <- as.matrix(table(out$Region,out$spp))

### mean number of ind per pop
#b <- aggregate(out$n,by=list(out$Region,out$spp),mean)
#colnames(b) <- c("Region","Spp","xbar_n")
md <- melt(out,id=c("Region","spp","pop"))
b <- cast(md,Region~spp,mean,na.rm=T)
b2 <- round(as.matrix(b),1)

b2[is.na(b2)] <- 0

pdf(paste("FINALFIGS/5_RandomForestCirclize/ALLSPECIES/Summary_sampleSize.pdf"),width=11,height=25)

tt = ttheme_default(colhead=list(fg_params=list(rot=90)))

table <- tableGrob(a,theme=tt)

title <- textGrob("Number of populations",gp=gpar(fontsize=20))
padding <- unit(2,"line")
table <- gtable_add_rows(table, 
                         heights = grobHeight(title) + padding,
                         pos = 0)
table <- gtable_add_grob(table, list(title),
                         t=1, l=1, 
                         r=ncol(table))
grid.draw(table)
grid.newpage()
table <- tableGrob(b2,theme=tt)

title <- textGrob("Mean# sample per population, averaged by region",gp=gpar(fontsize=20))
padding <- unit(2,"line")
table <- gtable_add_rows(table, 
                         heights = grobHeight(title) + padding,
                         pos = 0)
table <- gtable_add_grob(table, list(title),
                         t=1, l=1, 
                         r=ncol(table))
grid.draw(table)

dev.off()
