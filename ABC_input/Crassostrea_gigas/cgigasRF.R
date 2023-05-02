library(randomForest)
#library(caret)

metaname="cgiga_meta.csv"
species = "Cgigas"
meta = read.csv(metaname)
tab = read.table("SNPs_noZeros.txt",header=T)
rownames(tab)=tab[,1]
tab=tab[,-1]
gigas_newick=read.tree("~/GoogleDrive/data/Oyster/introduced/newickTrees/gigas_snp.tr")

pc = prcomp(tab)
tab.pred=predict(pc)

###comment out the following to use pca transform on data
tab.predict = tab
###

pop=sapply(strsplit(rownames(tab),"_"),function(x)x[1])

native_data = tab.pred[pop %in% meta$pop[meta$vars.reg=="1_Asia"],]
native_pops = as.factor(pop[pop %in% meta$pop[meta$vars.reg=="1_Asia"]])
intro_data =  tab.pred[!(pop %in% meta$pop[meta$vars.reg=="1_Asia"]),]
intro_pops =  as.factor(pop[!(pop %in% meta$pop[meta$vars.reg=="1_Asia"])])


rf = randomForest(x=native_data,y=native_pops)

pdf("gigas_random_forest.pdf")
plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()
