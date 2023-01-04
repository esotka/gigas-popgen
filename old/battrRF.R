library(strataG)
library(randomForest)
library(adegenet)
library(ape)
library(apex)
#library(caret)

metaname="battr_meta.csv"
species = "Gverm"
meta = read.csv(metaname)
gt = readRDS("battr_mtDNA_gtype")

##############haplo version
if (TRUE) {
aln=multidna2alignment(gt@sequences)
uhaplo=unique(aln$seq)
names(uhaplo) = paste0("hap",1:length(uhaplo))
dst = as.matrix(dist.dna(as.DNAbin(strsplit(uhaplo,""))))
hid = unname(sapply(aln$seq,function(x){names(uhaplo)[which(uhaplo==x)]}))
hapmat=matrix(0,nrow=aln$nb,ncol=length(uhaplo))
cnt=1
for(h in hid) 
{
    hapmat[cnt,as.numeric(gsub("hap","",h))]=1
    cnt=cnt+1
}
dmat = matrix(NA,nrow=nrow(hapmat),ncol=ncol(hapmat))
for (i in 1:nrow(hapmat))
{
    h = which(hapmat[i,]>0)
    dmat[i,]=dst[,h]    
}
tab=cbind(hapmat,dmat)

####################
} else {
########################make the data look like snps
gin = multidna2genind(gt@sequences)
tab = gin@tab
#######################
}

pop = getStrata(gt)

native_data = tab[pop %in% meta$source[meta$vars.reg=="1_Asia"],]
native_pops = as.factor(pop[pop %in% meta$source[meta$vars.reg=="1_Asia"]])
intro_data =  tab[!(pop %in% meta$source[meta$vars.reg=="1_Asia"]),]
intro_pops =  as.factor(pop[!(pop %in% meta$source[meta$vars.reg=="1_Asia"])])


rf = randomForest(x=native_data,y=native_pops)

pdf("battr_random_forest.pdf")
plt=data.frame(pred=predict(rf,newdata=intro_data),pop=intro_pops)
tbl=with(plt,table(pred,pop))
par(mar=c(7,4,4,4)+0.1)
heatmap(tbl,cexCol=1.7,cexRow=1.7)
heatmap(tbl[rowSums(tbl)>0,],cexCol=1.7,cexRow=1.7)
dev.off()
