library(strataG)
library(ape)
#library(randomForest)
library(ranger)
library(caret)

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

native = data.frame(cbind(native_pops,native_data))
introduced = data.frame(cbind(intro_pops,intro_data))
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
