##' take snp data from bcftool calls (vcf format, gzipped)
##' @param geno genotype structure created by modifying vcfR output
##' @return a rmetasim landscape with the same number of inds and genotypes as empirical data
##' @author allan strand
##' @export
create_gigas = function(geno)
{
    
    pops=sapply(strsplit(rownames(geno),"_"),function(x){x[1]})
    h=length(unique(pops))
    gl = list("0"=c(1,1),"1"=c(1,2),"2"=c(2,2))
        
    l =  landscape.new.empty()
    
    l <- landscape.new.intparam(l, h = h, s = 2)
    l <- landscape.new.switchparam(l, mp = 0)
    l <- landscape.new.floatparam(l)
    S <- matrix(c(0, 0, 1, 1), byrow = TRUE, nrow = 2)
    R <- matrix(c(0, 10, 0, 0), byrow = TRUE, nrow = 2)
    M <- matrix(c(0, 0, 0, 1), byrow = TRUE, nrow = 2)
    l <- landscape.new.local.demo(l, S, R, M)
    M <- R <- S <- matrix(0, nrow = h*2, ncol=h*2)
    
    l <- landscape.new.epoch(l, S = S, R = R, M = M, 
                             carry = rep(1000,h))
    
    imat = matrix(0,ncol=6,nrow=nrow(geno))
    imat[,4]=1:nrow(imat)
    
    
    for (i in 1:ncol(geno))
    {
        imat=cbind(imat,do.call(rbind,gl[as.character(geno[,i])])[,1:2])
        l <- landscape.new.locus(l, type = 2, ploidy = 2, 
                                 mutationrate = 0.000001, transmission = 0, numalleles = 2, 
                                 allelesize = 1, states=sample(c("A","T","C","G"),2,replace=F))
    }

    hm = matrix(0:((2*length(unique(pops))-1)),nrow=2)
    popids=hm[2,]
    names(popids)=unique(pops)
    imat[,1] <- popids[pops]
    l$individuals = imat
    l$intparam$nextid = max(imat[,3])+1
    is.landscape(l)
    l
}
