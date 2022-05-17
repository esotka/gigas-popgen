### PCA - downsampled Li + our Samples

rm(list=ls())
library(RColorBrewer)
library(scales) ## transparency with alpha
meta <- read.csv('data/MetaData_Cgigas.csv')
ind <- readLines('data/downsampled.713inds')
pop <- substr(ind,1,3)
reg <- as.character(meta$reg.code[match(pop,meta$PopID)])
reg[pop=="Chi"] <- "China"
reg[pop=="Jap"] <- "Jap"
reg[pop=="Kor"] <- "Kor"
reg[pop=="Sou"] <- "SouthAfrica"
#reg[substr(ind,1,10)=="China:Ning"] <- "Ningde"
reg <- substr(reg,1,3)
reg <- factor(reg)
dat <- read.delim('data/downsampled.ourSamples+Li.cov',header=F)
e <- eigen(dat)
pca1 <- e$vectors[,1]; pca2 <- e$vectors[,2]

pdf("output/PCA-downsampledLi+ourSamples.pdf")
plot(pca2~pca1,pch=20,col=alpha(brewer.pal(6,"Spectral")[reg],.5),xlim=range(pca1),ylim=range(pca2))
legend(x=0.0,y=-0.08,fill=brewer.pal(6,"Spectral"),legend=levels(reg),cex=.75)
dev.off()
