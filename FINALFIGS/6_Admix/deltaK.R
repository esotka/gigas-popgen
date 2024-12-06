
pdf("deltaK.pdf")
lnl = read.csv("admix_likelihood_reps.csv",header=F)
colnames(lnl) = c("rep","k","lnl")
lnl$k = as.numeric(gsub("k","",lnl$k))
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(x=lnl$k,y=lnl$lnl,xaxt="n",xlab="k",ylab="lnl")
mtext(at=2:8,2:8,side=1,line=1)
xbar <- tapply(lnl$lnl,lnl$k,mean)
std <- tapply(lnl$lnl,lnl$k,sd)
out <- data.frame(xbar,std)
# Evanno et al. 2005 ∆K = m(|L(K + 1) − 2 L(K ) + L(K − 1)|)/s[L(K )]
out$L.prime.k <- c(NA,xbar[-1]-xbar[-length(xbar)])
out$L.dblprime.k <- c(NA,out$L.prime.k[-c(1,length(xbar))]-(out$L.prime.k)[-(1:2)],NA)
out$delta <- out$L.dblprime.k/out$std
print(out)
plot(out$L.dblprime.k,xaxt="n",xlab="k",ylab="L.dblprime.k", type="b")
mtext(at=1:7,2:8,side=1,line=1)
#text(2,max(out$delta,na.rm=T)*.9,"C",cex=2) deltaK analysis
dev.off()
