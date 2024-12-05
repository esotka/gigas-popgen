### do crossvalidation
rm(list=ls())
library(caret)
library(ggplot2)
library(ranger)
load("gvermSNPforRangerCrossValidation.rda")

fitControl=trainControl(method="repeatedcv", number=10,repeats=10) #10-fold cv repeated 10 times
rangerGrid = expand.grid(mtry = round(seq(.1*ncol(native),0.95*ncol(native),length=5)),splitrule=c("gini","extratrees"),min.node.size=c(1,3,5,10) )
control = list(ranger=list(method="ranger",tuneGrid=rangerGrid))

rfits = lapply(control,function(x)
{
  print(x$method)
  train(reg~., data=native, method=x$method,
        tuneGrid=x$tuneGrid, trControl=fitControl)
})

save(rfits,file="gvermSNP_Ranger_out.rda")
#rf = ranger(reg~.,data=native, mtry=rfits[[1]]$bestTune$mtry,
#           splitrule=rfits[[1]]$bestTune$splitrule,
#            min.node.size=rfits[[1]]$bestTune$min.node.size)
#save(rf,file="gvermSNP_Ranger_out_rf.rda")

