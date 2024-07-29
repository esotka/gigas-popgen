library(strataG)
library(parallel)
library(ggplot2)
library(dplyr)

cores=12

path="osgresults/sep5/"

files = list.files(path=path, pattern="*.rda")
dat=data.frame(filename=files)
files=gsub("structout_","",gsub(".rda","",files))
df = data.frame(do.call(rbind,lapply(strsplit(files,"_"),function(x) x)))
names(df)=c("sourceA","sourceB","k","num.parents","gen.feral","burnin","chainlen","rep")
for (i in 3:ncol(df)) df[,i]=as.numeric(df[,i])
dat=cbind(df,dat)

if (!file.exists("rl.rda"))
    {
        rl = mclapply(1:length(files),mc.cores=cores,function(i) 
        {
            print(i)
            load(paste0(path,dat$filename[i]))
            overallFst=outlist$NE_vs_all_Fst[2,1]
            pwFst=outlist$pwFst
            pwFst[upper.tri(pwFst)] = 0
            pwFst = pwFst + t(pwFst)
            sp=outlist$suspectPops
            suspDiff=mean(pwFst[rownames(pwFst)%in%sp,!colnames(pwFst)%in%sp],na.rm=T)
            inSusp=mean(pwFst[rownames(pwFst)%in%sp,colnames(pwFst)%in%sp],na.rm=T)
            inOther=mean(pwFst[!rownames(pwFst)%in%sp,!colnames(pwFst)%in%sp],na.rm=T)
            het = data.frame(heterozygosity(outlist$g,by.strata=T,type="expected") %>%
                                    group_by(stratum) %>% summarise(He=mean(exptd.het)))
            rownames(het) = c("Normal","Suspect")
            
            res = data.frame(overallFst=overallFst,suspDiff=suspDiff,inSusp=inSusp,inOther=inOther,hetN=het[1,2],hetS=het[2,2])
            list(pwFst=pwFst,res=res,parms=as.data.frame(dat[i,]))
        })
        rl = rl[!sapply(rl,function(x) {class(x)%in%"try-error"})]        
        save(file="rl.rda",rl)
    } else {
        load("rl.rda")
    }
overall = do.call(rbind,lapply(rl,function(x){cbind(x$res,x$parms)}))
overall$pairing=paste(overall$sourceA,overall$sourceB,sep="_")

mns=overall %>% filter(sourceA=="CTRL") %>%mutate(diffRatio=inSusp/inOther)%>% select(overallFst,suspDiff,diffRatio) %>% colMeans()

png("among_pop_diversity.png")
### treating suspect and everything else as two populations (ignore hier)
overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, sourceA, sourceB, num.parents, gen.feral) %>% summarise(overallFst=mean(overallFst)) %>%
    mutate(sameSource=sourceA==sourceB) %>%
    ggplot(aes(y=overallFst,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    geom_hline(yintercept=mns[1],linetype="dashed", color = "red") + facet_wrap(~gen.feral) +
    ggtitle("Overall Fst among suspect and 'normal' pops\nTreating as if 2 populations")

dev.off()

pdf("among_pop_diversity.pdf")

### treating suspect and everything else as two populations (ignore hier)
overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, sourceA, sourceB, num.parents, gen.feral) %>% summarise(overallFst=mean(overallFst)) %>%
    mutate(sameSource=sourceA==sourceB) %>%
    ggplot(aes(y=overallFst,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    geom_hline(yintercept=mns[1],linetype="dashed", color = "red") + facet_wrap(~gen.feral) +
    ggtitle("Overall Fst among suspect and 'normal' pops\nTreating as if 2 populations")

overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, sourceA, sourceB, num.parents, gen.feral) %>% summarise(overallFst=mean(overallFst)) %>%
    mutate(sameSource=sourceA==sourceB) %>%
    ggplot(aes(y=overallFst,x=num.parents,group=pairing,color=sameSource)) + geom_point() + geom_line() +
    geom_hline(yintercept=mns[1],linetype="dashed", color = "red") + facet_wrap(~gen.feral) +
    ggtitle("Overall Fst among suspect and 'normal' pops\nTreating as if 2 populations")

###average pairwise between suspect and other populations
overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, num.parents, gen.feral) %>% summarise(suspDiff=mean(suspDiff)) %>%
    ggplot(aes(y=suspDiff,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    geom_hline(yintercept=mns[2],linetype="dashed", color = "red") + facet_wrap(~gen.feral) +
     ggtitle("Average pairwise difference between suspect pops and 'normal' pops")

overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, num.parents, gen.feral) %>% summarise(Diversity.Ratio=mean(inSusp/inOther)) %>%    
    ggplot(aes(y=Diversity.Ratio,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    geom_hline(yintercept=mns[3],linetype="dashed", color = "red") + facet_wrap(~gen.feral) +
    ggtitle("Ratio of average pairwise among suspect pops/average pairwise\namong 'normal' pops")


overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, num.parents, gen.feral) %>% summarise(He=mean(hetN)) %>%    
    ggplot(aes(y=He,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line()  + facet_wrap(~gen.feral) +
    ggtitle("Expected heterozygostiy in 'normal' pops")

overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, num.parents, gen.feral) %>% summarise(He=mean(hetS)) %>%    
    ggplot(aes(y=He,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    facet_wrap(~gen.feral) +
    ggtitle("Expected heterozygostiy in 'suspect' pops")

overall %>% filter(sourceA!="CTRL") %>% group_by(pairing, num.parents, gen.feral) %>% summarise(He=mean(hetS/hetN)) %>%    
    ggplot(aes(y=He,x=num.parents,group=pairing,color=pairing)) + geom_point() + geom_line() +
    facet_wrap(~gen.feral) +
    ggtitle("Expected heterozygosity ratio suspect/normal")


dev.off()

## structure plots
