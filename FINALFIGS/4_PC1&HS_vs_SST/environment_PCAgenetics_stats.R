### PCA environment vs PCA genetics
library(GGally)
pc_results <- read.csv("FINALFIGS/3_PCA/PC.out.csv")
meta <- read.csv("FINALFIGS/gigas_meta_41pop_env_FINAL.csv")
out <- data.frame(pc1=pc_results$PC1,meta[match(pc_results$Group.1,meta$pop),c(13:25)])
rownames(out) <- pc_results$Group.1
natnon <- meta$NatNon[match(pc_results$Group.1,meta$pop)]
## Correlations among variables with PCs
### Native and Introduced combined
### PC1 best correlates with 
#1) SSTmean (.773)  
#2) Mean Air Temp (.693)  
#3) SSTmin (0.681)  

pdf("FINALFIGS/4_PC1&HS_vs_SST/environment_PCAgenetics_stats.pdf",height=10,width=10)
print(ggscatmat(data=out,alpha=.1)) # PC1 and biotic  variables
dev.off()
print(sort(cor(out)[1,]))
heatmap(abs(cor(out)),main="native and introduced")



### Native only 
### PC1 best correlates with 
#1) Max air temp (0.787)
#2) SST mean (0.777)
#3) Min Air Temp (0.751)


print(sort(cor(out[natnon=="Native",])[1,]))
heatmap(abs(cor(out[natnon=="Native",])),main="native only")


### Introduced only 
### PC1 best correlates with 
#1) SST mean (0.744)
#2) SST min (0.714)
#3) Min air temp (0.707)


print(sort(cor(out[natnon=="Introduced",])[1,]))
heatmap(abs(cor(out[natnon=="Introduced",])),main="introduced only")

## plot of SSTmean (seen in all three comparisons) that best correlate with PC1

mysettings <- trellis.par.get()
mysettings$strip.background$col <- "lightgrey"
mysettings$box.dot$cex <- 0
trellis.par.set(mysettings)
xyplot(pc1~BO_sstmean | natnon,data=out,type=c("p","r"),main="r=0.773 (all); 0.777 (native); 0.744 (introduced)")
print(cor.test(~pc1+BO_sstmean,data=out))
print(cor.test(~pc1+BO_sstmean,data=out[natnon=="Native",]))
print(cor.test(~pc1+BO_sstmean,data=out[natnon=="Introduced",]))


