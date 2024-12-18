###  The idea is to get the RFCompAll files from the individual species data folders and then
###  parse it and compare the Bayes factors you get from Admix vs Shipping
###
###
library(tidyr)
library(dplyr)
library(ggplot2)

#cmd=paste0('find . -name "RFCompAllReps*" -print|grep priors')
#files=system(cmd,intern=TRUE)
#ret = do.call(rbind,lapply(files,function(fn) read.csv(fn)))
#write.table(file="RFCompAllReps.csv",sep=",",row.names=F,ret)

ret = read.csv("ABC_input/RFCompAllReps.csv")
ret$priorcat=basename(ret$priorfile)
ret$postship=ifelse(ret$chosen=="g2",
         ifelse(ret$comparison!="AdmixVML",1-ret$post,NA),
         ifelse(ret$comparison!="AdmixVML",ret$post,NA))
ret$postoyster=1-ret$postship
ret=ret[complete.cases(ret),]

##now there is a dataframe with all species and all comparisons, let's format for
##analysis

#head(ret)
wide = ret %>% filter(type=="untrans") %>%
    select(comparison,species,priorcat,priorShip,postship,postoyster) %>%
    pivot_wider(id_cols=c(species,priorcat),values_from=c(5:6),names_from=c(1))%>%
    mutate(dev=abs(postship_ShipVadmix-postship_ShipVML)) %>%
    arrange(priorcat,-dev)

head(data.frame(wide))

pdf("ABC_input/ML_versus_Admix_shipping_prob.pdf")
ggplot(wide,aes(x=postship_ShipVML,y=postship_ShipVadmix)) +
    geom_point() +
    geom_abline(slope=1)+facet_wrap(~priorcat)+
    xlab("Random Forest (ML)-based oyster model") +
    ylab("Admixture-based oyster model")+
    ggtitle("posterior probability of shipping")
dev.off()


png("ABC_input/ML_versus_Admix_shipping_prob_FLAT.png",height = 6,width = 6,units="in",res=400)
flat = wide[wide$priorcat=="flat_priors.csv",]
flat$sppAbb = c("Up","Mc","Af","Ba","Gv","Dv","HL6","HL1","Ht","Up2","Hj","Ph","Pm","Hs")
plot(x=flat$postship_ShipVML,y=flat$postship_ShipVadmix, xlab = "Random Forest (ML)-based oyster model", 
    ylab = "Admixture-based oyster model", 
    main = "posterior probability of shipping",cex=4)
segments(-1,-1,1,1,lty="dotted")
text(x=flat$postship_ShipVML,y=flat$postship_ShipVadmix,flat$sppAbb)
dev.off()

print(cor.test(x=flat$postship_ShipVML,y=flat$postship_ShipVadmix))
flat2 = flat[!flat$sppAbb=="Up",]
print(cor.test(x=flat2$postship_ShipVML,y=flat2$postship_ShipVadmix))
