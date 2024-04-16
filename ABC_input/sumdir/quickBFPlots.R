###
### quick visualization of BF plots
###
library(ggplot2);library(tidyr); library(dplyr)

###dat = read.csv("rfBayesFactor2.csv")
dat=do.call(rbind,lapply(list.files(path=".",pattern="rfBayesFactor.*csv"),read.csv))
dat$species=gsub("-","",dat$species)

pdf("distribution_of_priors.pdf")

ggplot(dat[dat$type=="untrans",],aes(x=species,y=(1-plotpost))) + geom_boxplot()+geom_jitter(width=0.2,size=0.3) +facet_grid(rows=vars(prior),cols=vars(type))+ylab("Posterior probability of shipping")+xlab("")+ggtitle("Distributions of posteriors")+theme(axis.text.x=element_text(angle=67,vjust=1,hjust=1,size=13))

dev.off()

dirSpecies = read.csv("dir_species_map.csv")
dirSpecies$species=gsub(" ","",dirSpecies$species)

prior = do.call(rbind,lapply(list.files(path=".",pattern="*prior.*csv"),function(fn)
{
    p=read.csv(fn)
    p$prior=gsub(".csv","",fn)
    p
}))
prior$priorShip = prior$probShip
prior = prior[,names(prior)!="probShip"]

prior = merge(prior,dirSpecies)

dat = dat%>%
#    filter(!grepl("SNP",species),!grepl("SSR",species))%>%
    group_by(species,type,prior) %>%
    summarise(oob=mean(oob),plotpost=mean(plotpost),lower.quartile=quantile(plotpost,0.25),upper.quartile=quantile(plotpost,0.75)) %>%
    mutate(postOyster=plotpost,postShip=(1-plotpost))

dat = merge(dat,prior)

long = pivot_longer(dat,cols=c("postOyster","postShip","priorShip")) %>% arrange(species,prior)
oy = long %>% filter(name=="postOyster",prior=="flat_priors",type=="untrans") %>% group_by(species) %>% summarize(mnPost=mean(value))
long = merge(oy,long)

long$species = reorder(as.factor(long$species),long$mnPost)

long=long[long$type=="untrans",]

pdf("compare_posteriors_and_priors.pdf")
ggplot(long[long$name!="priorShip",],aes(x=species,y=value,group=name,fill=name)) + geom_bar(stat="identity",position="dodge")+facet_grid(rows=vars(prior),cols=vars(type))+theme(axis.text.x=element_text(angle=67,vjust=1,hjust=1,size=13))+ylab("Posterior probability")+xlab("")+ggtitle("Compare the posteriors for Oyster and Shipping")+scale_fill_manual(values=c("red","blue"))

ggplot(long[long$name!="postOyster",],aes(x=species,y=value,group=name,fill=name)) + geom_bar(stat="identity",position="dodge")+facet_grid(rows=vars(prior))+theme(axis.text.x=element_text(angle=67,vjust=1,hjust=1,size=13))+ylab("Posterior probability")+xlab("")+ggtitle("Compare the prior for shipping to the posterior across species")+scale_fill_manual(values=c("blue","orange"))
dev.off()
