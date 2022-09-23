### PCA environment vs Hs
rm(list=ls())
library(GGally)
stats <- read.csv("FINALFIGS/2_basicStats/BasicStats_all.csv")
stats <- stats[complete.cases(stats$Hs),] # remove YOJ (n=1)
stats <- stats[!stats$n < 5,] # remove MAI (n = 4)
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
out <- data.frame(hs=stats$Hs,meta[match(stats$X,meta$pop),c(13:25)])
out$lat <- abs(meta$Latitude[match(stats$X,meta$pop)])
rownames(out) <- stats$X
natnon <- meta$NatNon[match(stats$X,meta$pop)]
## Correlations among variables with PCs
### Native and Introduced combined
### PC1 best correlates with 
#1) SSTmax (.549)  
#2) SSTmean (.535)  
#3) Max air temp (0.489)  
# lat is 5th

pdf("FINALFIGS/4_PC1&HS_vs_SST/environment_Hs_stats.pdf",height=10,width=10)
print(ggscatmat(data=out,alpha=.1)) # PC1 and biotic  variables
dev.off()
print(sort(abs(cor(out))[1,]))
heatmap(abs(cor(out)),main="native and introduced")



### Native only 
### PC1 best correlates with 
#1) lat (0.740)
#2) SST min (0.728)
#2) SST mean (0.725)
#4) Mean Air Temp (0.710)

print(sort(abs(cor(out[natnon=="Native",]))[1,]))
heatmap(abs(cor(out[natnon=="Native",])),main="native only")


### Introduced only 
### PC1 best correlates with 
#1) SST max (0.396)
#2) SST mean (0.343)
#3) Mean air temp (0.305)

print(sort(abs(cor(out[natnon=="Introduced",]))[1,]))
heatmap(abs(cor(out[natnon=="Introduced",])),main="Introduced only")

## plot of SSTmean (seen in all three comparisons) that best correlate with PC1

mysettings <- trellis.par.get()
mysettings$strip.background$col <- "lightgrey"
mysettings$box.dot$cex <- 0
trellis.par.set(mysettings)
print(xyplot(hs~BO_sstmean | natnon,data=out,type=c("p","r"),main="r=0.535 (all); 0.725 (native); 0.343 (introduced)"))
print(cor.test(~hs+BO_sstmean,data=out))
print(cor.test(~hs+BO_sstmean,data=out[natnon=="Native",]))
print(cor.test(~hs+BO_sstmean,data=out[natnon=="Introduced",]))


