library(strataG)
library(ape)
library(randomForest)
#library(caret)

metaname="gverm_meta.csv"
gname = "gverm_usat_gtype"
species = "Gverm"
meta = read.csv(metaname)
gin = readRDS(file=gname)
gigas_newick=read.tree("~/GoogleDrive/data/Oyster/introduced/newickTrees/gigas_snp.tr")

##get rid of individuals with NA strata,

#gin = gin[!(is.na(getStrata(gin))),,]

pop=sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
srcdf = data.frame(gigas_source=c("miy","miy","hon","hok","tok","tok","sea"),
                   reg=c("EU","PNW","hon","hok","tok","Cali","sea"))
gsdf=merge(data.frame(reg=getStrata(gin),origorder=1:length(getStrata(gin))),srcdf,all=T)
gsdf=gsdf[order(gsdf$origorder),]
sch = data.frame(reg=getStrata(gin),
                 intro=grepl("Introduced",getIndNames(gin)),
                 pop=pop,
                 ind=sapply(strsplit(getIndNames(gin),"_"),function(x) x[4]),
                 gigas_source = gsdf$gigas_source
                 )


############################


data_mat = gtypes2genind(gin)@tab

native_data = as.data.frame(data_mat[!sch$intro,])
native_pops = as.factor(sch$pop[!sch$intro])
intro_data =  as.data.frame(data_mat[sch$intro,])
intro_pops =  as.factor(sch$pop[sch$intro])


rf = randomForest(x=native_data,y=native_pops)

plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))

pdf("gverm_random_forest.pdf")
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()
