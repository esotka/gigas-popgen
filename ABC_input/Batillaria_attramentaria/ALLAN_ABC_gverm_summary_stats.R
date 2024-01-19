###
### summary stat calculator
###



summary_stats.new = function(gin,meta)
{

    meta.orig=meta
    rownames(meta)=meta$longpop
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"

    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
    } else {
        longpops=sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
    }
    

names(longpops) = getIndNames(gin)
orig.strata = getStrata(gin)
setStrata(gin) = longpops
popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]

names(popHet) = paste0(meta$pop,"Het")

popOverall = overallTest(gin,nrep=0)$result[,1]
names(popOverall) = paste0(names(popOverall),".pop")
#popPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
#popPair=sapply(popPW,function(x) x$result["Fst",1])
#names(popPair) = sapply(popPW,function(x) paste(names(x$strata.freq),collapse="_"))
    
intro = meta[longpops,"intro"]
names(intro) = getIndNames(gin)
setStrata(gin) = intro
intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))

reg=meta[longpops,"gigas_source"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(regHet) =  unique(reg)
regM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(regM) =  unique(reg)
regOverall = overallTest(gin,nrep=0)$result[,1]
names(regOverall) = paste0(names(regOverall),".reg")
regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
regPair=sapply(regPW,function(x) x$result["Fst",1])
names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))



reg=meta[longpops,"reg"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
reg2Het = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
reg2M = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
reg2Overall = overallTest(gin,nrep=0)$result[,1]
names(reg2Overall) = paste0(names(reg2Overall),".reg2")
reg2PW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
reg2Pair=sapply(reg2PW,function(x) x$result["Fst",1])
names(reg2Pair) = sapply(reg2PW,function(x) paste(names(x$strata.freq),collapse="_"))


c(overallHet,popHet,intro_vs_nativeHet,regHet,reg2Het,
    popOverall,intro_vs_nativeOverall,regOverall,reg2Overall,
  intro_vs_nativePair,regPair,reg2Pair
  )
}


summary_stats.orig = function(gin,meta)
{

    meta.orig=meta
    rownames(meta)=meta$longpop
    overallHet = mean(heterozygosity(gin)[,2])
    names(overallHet)="overallHet"

    if (sum(grepl("/",getIndNames(gin)))==length(getIndNames(gin)))
    {
        longpops=getStrata(gin)
    } else {
        longpops=sapply(strsplit(getIndNames(gin),"_"),function(x) x[3])
    }
    

names(longpops) = getIndNames(gin)
orig.strata = getStrata(gin)
setStrata(gin) = longpops
popHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
popM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(popHet) = paste0(meta$pop,"Het")
names(popM) = paste0(meta$pop,"M")
popOverall = overallTest(gin,nrep=0)$result[,1]
names(popOverall) = paste0(names(popOverall),".pop")
#popPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
#popPair=sapply(popPW,function(x) x$result["Fst",1])
#names(popPair) = sapply(popPW,function(x) paste(names(x$strata.freq),collapse="_"))
    
intro = meta[longpops,"intro"]
names(intro) = getIndNames(gin)
setStrata(gin) = intro
intro_vs_nativeHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(intro_vs_nativeHet) = c("nativeHet","introducedHet")
intro_vs_nativeM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(intro_vs_nativeHet) = c("nativeM","introducedM")
intro_vs_nativeOverall = overallTest(gin,nrep=0)$result[,1]
names(intro_vs_nativeOverall) = paste0(names(intro_vs_nativeOverall),".intro")
introPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
intro_vs_nativePair=sapply(introPW,function(x) x$result["Fst",1])
names(intro_vs_nativePair) = sapply(introPW,function(x) paste(names(x$strata.freq),collapse="_"))

reg=meta[longpops,"gigas_source"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
regHet = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
names(regHet) =  unique(reg)
regM = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
names(regM) =  unique(reg)
regOverall = overallTest(gin,nrep=0)$result[,1]
names(regOverall) = paste0(names(regOverall),".reg")
regPW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
regPair=sapply(regPW,function(x) x$result["Fst",1])
names(regPair) = sapply(regPW,function(x) paste(names(x$strata.freq),collapse="_"))



reg=meta[longpops,"reg"]
names(reg) = getIndNames(gin)
setStrata(gin) = reg
reg2Het = with(heterozygosity(gin,by.strata=T),aggregate(cbind(exptd.het=exptd.het),by=list(stratum=stratum),mean))[,2]
reg2M = with(mRatio(gin,by.strata=T),aggregate(cbind(m.ratio=m.ratio),by=list(stratum=stratum),mean,na.rm=T))[,2]
reg2Overall = overallTest(gin,nrep=0)$result[,1]
names(reg2Overall) = paste0(names(reg2Overall),".reg2")
reg2PW=pairwiseTest(gin,nrep=0,stats="Fst_prime",quiet=T)
reg2Pair=sapply(reg2PW,function(x) x$result["Fst",1])
names(reg2Pair) = sapply(reg2PW,function(x) paste(names(x$strata.freq),collapse="_"))


c(overallHet,popHet,intro_vs_nativeHet,regHet,reg2Het,
  popM,intro_vs_nativeM,regM,reg2M,
  popOverall,intro_vs_nativeOverall,regOverall,reg2Overall,
  intro_vs_nativePair,regPair,reg2Pair
  )
}


######

summary_stats = summary_stats.new
