###
### quick visualization of BF plots
###
library(ggplot2);library(tidyr); library(dplyr)

dat = read.csv("rfBayesFactor.csv")
dat$species=gsub("-","",dat$species)
dat$postOyster=dat$plotpost
dat$postShip=1-dat$plotpost



long = pivot_longer(dat[,c(1,5,6,10,11)],cols=c("postOyster","postShip")) %>% arrange(species,prior)
oy = long %>% filter(name=="postOyster",prior=="flat_priors") %>% group_by(species) %>% summarize(mnPost=mean(value))
long = merge(oy,long)
long$species = reorder(as.factor(long$species),long$mnPost)

ggplot(long,aes(x=species,y=value,group=name,color=name)) + geom_bar(stat="identity",position="dodge")+facet_grid(rows=vars(prior))+theme(axis.text.x=element_text(angle=45,vjust=0.75,hjust=1))
