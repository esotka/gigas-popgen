library(strataG)
library(ape)
#library(randomForest)
library(ranger)
library(caret)

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
tab.pred = tab
###

pop=sapply(strsplit(rownames(tab),"_"),function(x)x[1])

native_data = tab.pred[pop %in% meta$pop[meta$vars.reg=="1_Asia"],]
native_pops = as.factor(pop[pop %in% meta$pop[meta$vars.reg=="1_Asia"]])
intro_data =  tab.pred[!(pop %in% meta$pop[meta$vars.reg=="1_Asia"]),]
intro_pops =  as.factor(pop[!(pop %in% meta$pop[meta$vars.reg=="1_Asia"])])

native = data.frame(cbind(native_pops,data.frame(native_data)))
introduced = data.frame(cbind(intro_pops,data.frame(intro_data)))
names(native)[1]="pop"


names(introduced)[1]="pop"

########################## training and testing

intrain = createDataPartition(native$pop,p=0.8,list=F)

train=native[intrain,]
test =native[-intrain,]

fitControl=trainControl(method="repeatedcv", number=10,repeats=10) #10-fold cv repeated 10 times
rangerGrid = expand.grid(mtry = round(seq(5,50,length=5)),splitrule=c("gini","extratrees"),min.node.size=c(1,3,5,10) )

control = list(ranger=list(method="ranger",tuneGrid=rangerGrid))
               
rfits = lapply(control,function(x)
{
    print(x$method)
        train(pop~., data=train, method=x$method,
               tuneGrid=x$tuneGrid, trControl=fitControl)
    })

rf = ranger(pop~.,data=native, mtry=rfits[[1]]$bestTune$mtry,
            splitrule=rfits[[1]]$bestTune$splitrule,
            min.node.size=rfits[[1]]$bestTune$min.node.size)

pred = predict(rf,data=introduced)$predictions
actual=introduced$pop

tbl = table(pred,actual)
