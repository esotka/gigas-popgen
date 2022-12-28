###
### R script to generate PCA from a covariance matrix stored on disk
### Uses a bamfile and metadata to classify points, this happens in setup_meta.R
###

source("setup_meta.R")

p="."
f=list.files(path=p,pattern="*.qopt")

#######################################################
plotadmix = function(m=as.matrix(adin[[1]]$q),
	 cldf=iddf,
	 cp=1:nrow(m),
	 sortcols=c("Region","pop"),
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
   barplot(pm,beside=F,border=NA,col=c(cp[1]:cp[ncol(m)]),space=0,horiz=T,axes=F,names.arg=rep("",length(plotdf$id)),xpd=NA,
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
      mtext(side=2,text=r1_coords$vals,at=fac*rowMeans(r1_coords[,c("strt","stp")]))
      for(j in 1:nrow(r1_coords))
      {
        x=-1*c(0.01,0.012)[(1+(j%%2))]
        segments(x,r1_coords[j,"strt"],x,r1_coords[j,"stp"],type="l", lwd=2,col="black")   
      }

      }	

      if (right)
      {
      fac = 1 #-inds * 1.5      
      mtext(side=4,text=r2_coords$vals,at=fac*rowMeans(r2_coords[,c("strt","stp")]),cex=0.7,adj=-0.4)
      for(j in 1:nrow(r2_coords))
      {
        x=c(1.003,1.008)[(1+(j%%2))]
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

pdf("admix.pdf",width=14,height=10)

cp=1:15
cldf=iddf

layout(matrix(c(1:length(adin),rep(length(adin)+1,length(adin))),ncol=length(adin),nrow=2,byrow=T),
	widths=c(2.5,rep(1,length(adin)-2),2),
	heights=c(10,3))
for (i in 1:length(adin))
{
    m=as.matrix(adin[[i]]$q)
    if (i==1)
        plotadmix(m,cldf,left=T,right=F) else if (i==length(adin)) plotadmix(m,cldf,left=F,right=T)
                                         else plotadmix(m,cldf,left=F,right=F) 

}

par(mar=c(5,8,1,2))
print(diff(sapply(adin,function(x) x$LL)))

plot(diff(sapply(adin,function(x) x$LL))~sapply(adin,function(x) x$k)[-1],
     xlim=c(1,length(adin)+1.68),
     type="b",lwd=2,cex=1.5,pch=16,axes=F,
     xlab="",ylab="")
axis(2)
par(las=3)
mtext(line=4,side=2,text="Log-Likelihood for admixture model",cex=0.7)
legend(x=1,y=200000,legend="Difference in LogLike among previous and current K")
dev.off()




