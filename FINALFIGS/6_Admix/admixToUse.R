###
### R script to generate PCA from a covariance matrix stored on disk
### Uses a bamfile and metadata to classify points, this happens in setup_meta.R
### plots individual assignments from the qopt files, but gets the likelihoods for each K
### in a dataframe
###
library(Cairo)
library(dplyr)
rm(list=ls())
source("setup_meta.R")

p="."      #"FINALFIGS/6_Admix/"
f=list.files(path=p,pattern="*.qopt")
likes = read.csv("admix_likelihood_reps.csv",header=F)
names(likes) = c("rep","k","llike")
likes$k = as.numeric(gsub("k","",likes$k))
sumlike=likes %>% arrange(rep,k) %>% group_by(rep) %>% reframe(rep=rep,d=c(NA,diff(llike))) %>%
    mutate(k=rep(c("","2-3","3-4","4-5","5-6","6-7","7-8"),50)) %>% filter(!is.na(d)) %>%
    group_by(k) %>% summarise(mnd=mean(d),mind=min(d),maxd=max(d),
                              sdd=sd(d),n=n(),se=sd(d)/sqrt(n()),
                              percent05=quantile(d,0.05),
                              percent95=quantile(d,0.95))
sumlike=as.data.frame(sumlike)
sumlike$plotk=3:8

kcol <- c(       "black", #1
                 "red",#2
                 "darkgreen",#3
                 "gainsboro",#4
                 "yellow",#5
                 "deepskyblue",#6
                 "brown",#7
                 "dodgerblue4")#8
colorder <- list(
                   c(1,2), #ks=2
                   c(2,1,3), #ks=3
                   c(3,4,1,2), #ks=4
                   c(2,5,4,3,1), #ks=5
                   c(5,2,6,1,3,4), #ks=6
                   c(6,2,7,1,3,4,5), #ks=7
                   c(7,2,1,4,3,6,5,8)) #ks=8


#######################################################
plotadmix = function(m=as.matrix(adin[[1]]$q),
	 cldf=iddf,
	 cp=1:nrow(m),
	 sortcols=c("Region2","pop"),
	 left=TRUE,
	 right=TRUE)	
{ 
    
   plotdf=data.frame(cbind(cldf,m))
   names(plotdf)[(ncol(cldf)+1):ncol(plotdf)] = paste0("q",1:ncol(m))
## Sorting
if ((!is.null(sortcols))&(length(sortcols)>0)) {
   if (length(sortcols)<2) rep(sortcols,2)
   plotdf = plotdf[order(plotdf[,sortcols[1]],-plotdf[,"Latitude"]),]
}  


par(las=2)
mar=c(0,0.1,1,0.1)+0.1
xlim=c(0,1)
if (left) {mar[2]=6 ; xlim[1]=-0.2}
if (right){ mar[4]=6; xlim[2]=1.02}

par(mar=mar)
   pm=t(as.matrix(plotdf[,(ncol(cldf)+1):ncol(plotdf)]))
   barplot(pm,beside=F,border=NA,col=c(kcol[colorder[[i]]]),space=0,horiz=T,axes=F,names.arg=rep("",length(plotdf$id)),xpd=NA,
           xlim=xlim,main=paste("K=",ncol(m)))
   
   if ((!is.null(sortcols))&(length(sortcols)>0)) {
      r1 = rle(as.character(plotdf[,sortcols[1]]))
      r2 = rle(as.character(plotdf[,sortcols[2]]))	

      r1_coords = data.frame(stp=cumsum(r1$lengths),
 			      vals=r1$values)
      r1_coords$strt = r1_coords$stp - r1$lengths

      r2_coords = data.frame(stp=cumsum(r2$lengths),
 			      vals=r2$values)
      r2_coords$strt = r2_coords$stp - r2$lengths

      
      inds=nrow(plotdf)

      if (left)
      {
      fac = 1 #-inds * 1.5      
      mtext(side=2,line=-2,text=r1_coords$vals,at=fac*rowMeans(r1_coords[,c("strt","stp")]))
      for(j in 1:nrow(r1_coords))
      {
        #x=-1*c(0.01,0.012)[(1+(j%%2))]
        x=-1*c(0.02,0.05)[(1+(j%%2))]
        segments(x,r1_coords[j,"strt"],x,r1_coords[j,"stp"],type="l", lwd=2,col="black")   
      }

      }	

      if (right)
      {
      fac = 1 #-inds * 1.5      
      mtext(side=4,line=.075,text=r2_coords$vals,at=fac*rowMeans(r2_coords[,c("strt","stp")]),cex=0.7,adj=-0.4)
      for(j in 1:nrow(r2_coords))
      {
        x=c(1.01,1.05)[(1+(j%%2))]
        segments(x,r2_coords[j,"strt"],x,r2_coords[j,"stp"],type="l", lwd=2,col="black")   
      }

      }	

   }  


}
########################################### end of plotting function declaration


##################### now use the function

adin <-  vector("list",length(f))
cnt=1
print(length(f))
for (x in f)
{
print(x)
    froot=gsub(".qopt","",x)
print(froot)
    k=as.numeric(gsub("admix_","",froot))
    l=readLines(paste0(p,"/",froot,".log"))
    l=gsub("best like=","",l[grep("best like",l)])
    LL=as.numeric(gsub(" after.*","",l))
    q=read.table(file=paste0(p,"/",x),header=F)
    adin[[cnt]] =  list(k=k,LL=LL,q=q)
    cnt=cnt+1
}

print("about to plot")

#pdf("admix.pdf",width=14,height=10)
CairoPNG("admix.png",width=1400,height=1000)

cp=1:15
cldf=iddf

layout(matrix(c(1:length(adin),rep(length(adin)+1,length(adin))),ncol=length(adin),nrow=2,byrow=T),
	widths=c(2.5,rep(1,length(adin)-2),2),
	heights=c(10,3))
for (i in 1:length(adin))
{
    m=as.matrix(adin[[i]]$q)
    if (i==1)
        plotadmix(m,cldf,left=T,right=F) else if (i==length(adin)) plotadmix(m,cldf,left=F,right=T) else plotadmix(m,cldf,left=F,right=F) 

}



par(mar=c(5,8,1,2))

plot(mnd~plotk,
     xlim=c(1,length(adin)+1.68),
     type="b",lwd=2,cex=2,pch=16,axes=F,
     ylim=c(min(percent05),max(percent95)),
     xlab="",ylab="",data=sumlike)
axis(2)
par(las=3)
mtext(line=4,side=2,text="delta Log-Like",cex=0.8)
for (i in 1:nrow(sumlike))
    with(sumlike,points(x=c(plotk[i],plotk[i]),y=c(percent05[i],percent95[i]),type="l"))
#legend(x=1,y=300000,legend="Difference in LogLike among previous and current K")
dev.off()




