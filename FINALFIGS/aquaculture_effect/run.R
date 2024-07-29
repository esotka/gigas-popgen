rm(list=ls())

cores=12
library(parallel)
library(rmetasim)
library(vcfR)
library(strataG)

source("filterVCFfunction.R")
source("create_landscape.R")

num.parents = 2  #number of parents chosen from each source population to found new ones
gen.after.release=10 #number of generations after release



geno=filterVCF(filelist=c("all.bamlist_NC_047559.1.vcf.gz"))
#geno=filterVCF(filelist= c("all.bamlist_NC_047559.1.vcf.gz", "all.bamlist_NC_047560.1.vcf.gz", "all.bamlist_NC_047561.1.vcf.gz", "all.bamlist_NC_047562.1.vcf.gz", "all.bamlist_NC_047563.1.vcf.gz", "all.bamlist_NC_047564.1.vcf.gz", "all.bamlist_NC_047565.1.vcf.gz", "all.bamlist_NC_047566.1.vcf.gz", "all.bamlist_NC_047567.1.vcf.gz", "all.bamlist_NC_047568.1.vcf.gz"))

pops = sapply(strsplit(rownames(geno),"_"),function(x) x[[1]][1])
popids=(1:length(unique(pops)))-1
names(popids)=unique(pops)
names(popids)[which(names(popids)=="QB2")]="QUI"


suspectPops=c("BBH","TRO","WAD","SBB","CHI","SVA","LOS")
sourceA="SAM"
sourceB="SAM"

l = create_gigas(geno)

l.mod = l.orig = l

ind=l.mod$individuals
ind=ind[which((landscape.populations(l)-1)%in%popids[c(sourceA,sourceB)]),]
ind2=NULL
for (cl in unique(ind[,1]))
{
    ind2=rbind(ind2,ind[ind[,1]==cl,][sample(1:length(which(ind[,1]==cl)),num.parents),])
}
ind=ind2


l.mod$individuals=ind

#set up the reproduce into the suspect pops
Rchng=as.matrix(expand.grid(D=(popids[suspectPops]*2+1),S=(popids[c(sourceA,sourceB)]*2+2)))
l.mod$demography$epochs[[1]]$R[as.matrix(Rchng)]=100  #can be anything for number of offspring

#set up the game transfer from source pops
Mchng=unique(as.matrix(expand.grid(S=(popids[c(sourceA,sourceB)]*2+2),D=(popids[c(sourceA,sourceB)]*2+2))))
Mchng=Mchng[Mchng[,1]!=Mchng[,2],]
l.mod$demography$epochs[[1]]$M[as.matrix(Mchng)]=1

##don't allow exchange of gametes within pops for the initial cross unless same popsource
l.mod$demography$localdem[[1]]$LocalM[2,2]=ifelse(sourceA==sourceB,1,0)
##don't allow repro within pops for the initial cross
l.mod$demography$localdem[[1]]$LocalR[1,2]=0

##simulate a generation of mixing between the sources to form the suspect pops
l.mod2=landscape.simulate(l.mod,1)

##return the demography to closer to original to prevent any more movement
l.mod2$demography$epochs[[1]]$R[as.matrix(Rchng)]=0  
l.mod2$demography$epochs[[1]]$M[as.matrix(Mchng)]=0

l.mod2$demography$localdem[[1]]$LocalM=matrix(c(0,0,0,1),ncol=2)
l.mod2$demography$localdem[[1]]$LocalR=matrix(c(0,0,5.4,0),ncol=2)
l.mod2$demography$localdem[[1]]$LocalS=matrix(c(0,0.09,0,0.6),ncol=2)

l.mod2=landscape.simulate(l.mod2,gen.after.release)
l.mod2$individuals[,3]=0 #make everybody have the same time as the rest of l

#now calculate the size of the original suspect populations so we can sample from l.mod2 and replace l
tbl=table(landscape.populations(l))
suspectSizes=tbl[which((as.numeric(names(tbl))-1) %in% popids[suspectPops])]
names(suspectSizes)=as.numeric(names(suspectSizes))-1
print(dim(l$individuals))
l.bak=l
l$intparam$curentgen=l.mod2$intparam$curentgen+1



for (pid in names(suspectSizes))
{
  print(suspectSizes[pid])
    print(as.numeric(pid)+1)

    l$individuals = l$individuals[landscape.populations(l)!=(as.numeric(pid)+1),]

    print(dim(l$individuals))
    ind=l.mod2$individuals[landscape.populations(l.mod2)==(as.numeric(pid)+1),] #individuals from this population, but admixed source

    print(dim(ind))
#    ind[,1]=ind[,1]+1 
    l$individuals=rbind(l$individuals,ind[sample(1:nrow(ind),suspectSizes[pid],F),]) #sample the same number as were collected
    print(dim(l$individuals))
}

l$individuals = l$individuals[order(l$individuals[,1]),]
print(dim(l$individuals))
table(landscape.populations(l))

l.admix=l

## should have a popids object with names equal to pops with plotting order prepended:
meta=read.csv("gigas_bam_qualities.csv")
names(popids)%in%meta$truepop
meta=meta[order(meta$Region,meta$Latitude),]
popidsort=unique(meta[meta$truepop%in%names(popids),c("truepop")])
popiddf=data.frame(name=popidsort,ord=1:length(popidsort))
popiddf=popiddf[order(popiddf$name),]
names(popids)=paste(sprintf("%02i",popiddf$ord),names(popids),sep="_")



##slightly change the plotting order to group the suspect and source populations together
names(popids)[grep(sourceA,names(popids))]=paste("97",sourceA,sep="_")
names(popids)[grep(sourceB,names(popids))]=paste("98",sourceB,sep="_")
for (sp in suspectPops) names(popids)[grep(sp,names(popids))]=paste("99",sp,sep="_")


g.admix=landscape2gtypes(l.admix)
g.orig=landscape2gtypes(l.orig)

krange=3:6
nreps=3
burnin=1000
numreps=1000


parms = expand.grid(k=krange,reps=1:nreps)
parms=parms[sample(1:nrow(parms),nrow(parms),F),] #hopefully helps with randomization

sr.orig=mclapply(1:nrow(parms),mc.cores=cores,mc.preschedule=F,
                 function(r)
                 {
                     Sys.sleep(runif(1,0,1))
                     k=parms[r,"k"]
                     res=structureRun(g.orig,k.range=k:k,num.k.rep=1,burnin=burnin,numreps=numreps,pops=names(popids),seed=floor(runif(1)),label=paste(parms[r,"k"],parms[r,"reps"],sep="_"))[[1]]
                     res$q.mat$orig.pop=names(popids)[as.numeric(res$q.mat$orig.pop)]
                     res

                 })
names(sr.orig)=paste(parms[,1],parms[,2],sep="_")
class(sr.orig)=c("structure.result", "list")



sr.admix=mclapply(1:nrow(parms),mc.cores=cores,mc.preschedule=T,
                  function(r)
                  {
                      Sys.sleep(runif(1,0,1))
                      k=parms[r,"k"]
                      res=structureRun(g.admix,k.range=k:k,num.k.rep=1,burnin=burnin,numreps=numreps,pops=names(popids),seed=floor(runif(1)),label=paste(parms[r,"k"],parms[r,"reps"],sep="_"))[[1]]
                      res$q.mat$orig.pop=names(popids)[as.numeric(res$q.mat$orig.pop)]
                      res
                  })
names(sr.admix)=paste(parms[,1],parms[,2],sep="_")
class(sr.admix)=c("structure.result", "list")


save(file="structout.rda",sr.admix,sr.orig,geno,l.admix,l.orig,g.admix,g.orig)


pdf (paste("structout_",num.parents,"_",gen.after.release,".pdf"),height=14,width=10)

evanno(sr.admix)
evanno(sr.orig)

q.admix=clumpp(sr.admix,k=4)
structurePlot(q.admix)

q.orig=clumpp(sr.orig,k=4)
structurePlot(q.orig,horiz=T)

dev.off()
