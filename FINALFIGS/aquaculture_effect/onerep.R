
args=commandArgs(trailingOnly=TRUE)
if (length(args)<1)
{
   r=1
} else r=as.numeric(args)+1

library(parallel)
library(rmetasim)
library(vcfR)
library(strataG)

source("filterVCFfunction.R")
source("create_landscape.R")

structure=F
onechrom=TRUE


krange=5  #3:7
nreps=5

parms = expand.grid(k=krange,reps=1:nreps,srcA=c("CTRL","OKA","SAM","OHK"),
                    srcB=c("SAM","OKA","OHK"),burnin=c(3000),
                    numreps=c(17000),num.parents=c(2,4,6,8,10,15),
                    gen.after.release=c(20,10,5,3,7,1))

ctr=which((parms$srcA=="CTRL") & (parms$srcB!="SAM"))
if (length(ctr)>0) parms=parms[-ctr,]
parms[,c("srcA","srcB")]=(apply(parms[,c("srcA","srcB")],2,sort))
parms=unique(parms)
rownames(parms)=1:nrow(parms)

print(paste("parms ",nrow(parms)," long"))

suspectPops=c("BBH","TRO","WAD","SBB","CHI","SVA","LOS")
sourceA=parms$srcA[r]
sourceB=parms$srcB[r]

if (onechrom==TRUE)
{
    geno=filterVCF(filelist=c("all.bamlist_NC_047559.1.vcf.gz"))
} else {
    geno=filterVCF(filelist= c("all.bamlist_NC_047559.1.vcf.gz", "all.bamlist_NC_047560.1.vcf.gz", "all.bamlist_NC_047561.1.vcf.gz", "all.bamlist_NC_047562.1.vcf.gz", "all.bamlist_NC_047563.1.vcf.gz", "all.bamlist_NC_047564.1.vcf.gz", "all.bamlist_NC_047565.1.vcf.gz", "all.bamlist_NC_047566.1.vcf.gz", "all.bamlist_NC_047567.1.vcf.gz", "all.bamlist_NC_047568.1.vcf.gz"))
}

pops = sapply(strsplit(rownames(geno),"_"),function(x) x[[1]][1])
popids=(1:length(unique(pops)))-1
names(popids)=unique(pops)
names(popids)[which(names(popids)=="QB2")]="QUI"

l = create_gigas(geno)

l.mod = l.orig = l

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


if (parms$srcA[r]=="CTRL") ##run on the original
{  
    print(paste("in CTRL, about to convert to gtpyes"))
    print(parms[r,])
    g.orig=landscape2gtypes(l.orig)

if (structure==T)
    k=parms[r,"k"]
    res=structureRun(g.orig,k.range=k:k,num.k.rep=1,burnin=parms$burnin[r],numreps=parms$numreps[r],pops=names(popids),seed=floor(runif(1)),label=paste(parms[r,"k"],parms[r,"reps"],sep="_"))[[1]]
    res$q.mat$orig.pop=names(popids)[as.numeric(res$q.mat$orig.pop)]
    sr.orig = res
    names(sr.orig)=paste(parms$k[r],parms$gen.after.release[r],sep="_")
    class(sr.orig)=c("structure.result", "list")
} else sr.orig=NULL

    sr=sr.orig
    g=g.orig
    l=l.orig


} else {  ##create and run on the admix

    ind=l.mod$individuals
    ind=ind[which((landscape.populations(l)-1)%in%popids[sapply(c(sourceA,sourceB),function(x) grep(x,names(popids)))]),]
    ind2=NULL
    for (cl in unique(ind[,1]))
    {
        ind2=rbind(ind2,ind[ind[,1]==cl,][sample(1:length(which(ind[,1]==cl)),parms$num.parents[r],replace=F),])
    }
    ind=ind2
    
    
    l.mod$individuals=ind
    
                                        #set up the reproduce into the suspect pops
    Rchng=as.matrix(expand.grid(D=(popids[sapply(suspectPops,function(x) grep(x,names(popids)))]*2+1),S=(popids[sapply(c(sourceA,sourceB),function(x) grep(x,names(popids)))]*2+2)))
    l.mod$demography$epochs[[1]]$R[as.matrix(Rchng)]=100  #can be anything for number of offspring
    
                                        #set up the game transfer from source pops
    Mchng=unique(as.matrix(expand.grid(S=(popids[sapply(c(sourceA,sourceB),function(x) grep(x,names(popids)))]*2+2),D=(popids[sapply(c(sourceA,sourceB),function(x) grep(x,names(popids)))]*2+2))))
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
    l.mod2$demography$localdem[[1]]$LocalR=matrix(c(0,0,7,0),ncol=2)
    l.mod2$demography$localdem[[1]]$LocalS=matrix(c(0,0.09,0,0.6),ncol=2)
    
    l.mod2=landscape.simulate(l.mod2,parms$gen.after.release[r])
    l.mod2$individuals[,3]=0 #make everybody have the same time as the rest of l
    
                                        #now calculate the size of the original suspect populations so we can sample from l.mod2 and replace l
    tbl=table(landscape.populations(l))
    suspectSizes=tbl[which((as.numeric(names(tbl))-1) %in% popids[sapply(suspectPops,function(x) grep(x,names(popids)))])]
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
        l$individuals=as.matrix(rbind(l$individuals,ind[sample(1:nrow(ind),suspectSizes[pid],replace=F),])) #sample the same number as were collected
        rownames(l$individuals)=1:nrow(l$individuals)
        print(dim(l$individuals))
    }
    
    l$individuals = l$individuals[order(l$individuals[,1]),]

    l$individuals[,4]=1:nrow(l$individuals)

    print(dim(l$individuals))
    table(landscape.populations(l))
    
    l.admix=l
    
    print(paste("in ADMIX, about to convert to gtypes"))
    print(parms[r,])

    
    g.admix=landscape2gtypes(l.admix)
    if (structure==T)
        {
    k=parms[r,"k"]
    res=structureRun(g.admix,k.range=k:k,num.k.rep=1,burnin=parms$burnin[r],numreps=parms$numreps[r],pops=names(popids),seed=floor(runif(1)),label=paste(parms[r,"k"],parms[r,"reps"],sep="_"))[[1]]
    res$q.mat$orig.pop=names(popids)[as.numeric(res$q.mat$orig.pop)]
    sr.admix=res
    
    names(sr.admix)=paste(parms$k[r],parms$gen.after.release[r],sep="_")
    class(sr.admix)=c("structure.result", "list")
    } else sr.admix=NULL
    sr=sr.admix
    g=g.admix
    l=l.admix
}

#calc some popgen stats
tmp = names(popids)[as.numeric(getStrata(g))]
tmp2 = sapply(strsplit(tmp,"_"),function(x) x[2])
names(tmp2)=getIndNames(g)
setStrata(g)=tmp2
pwFst=pairwiseMatrix(pairwiseTest(g,nrep=1),stat="Fst")

tmp3 = tmp2 %in% suspectPops
names(tmp3)=names(tmp2)
setStrata(g) = tmp3
NE_vs_all_Fst=pairwiseMatrix(pairwiseTest(g,nrep=1),stat="Fst")

outlist=list(sr=sr,g=g,k=parms$k[r],
		srcA=parms$srcA[r],srcB=parms$srcB[r],
		suspectPops=suspectPops,
             	numparents=parms$num.parents[r],
             	gens.after.release=parms$gen.after.release[r],
	     	pwFst=pwFst,NE_vs_all_Fst=NE_vs_all_Fst)

print(paste("r is equal to :",r))

save(file=paste0("structout_",parms$srcA[r],"_",ifelse(parms$srcA[r]=="CTRL","CTRL",parms$srcB[r]),"_",
                 parms$k[r],"_",parms$num.parents[r],"_",parms$gen.after.release[r],"_",parms$burnin[r],"_",parms$numreps[r],"_",
                 parms$rep[r],".rda"),
     outlist)
