landscape.sample.pops <- function(Rland,ns=NULL,pvec=NULL)
    {
        if (is.null(pvec))
            {
                Rland
            } else { #there are pops
                pops <- 1:Rland$intparam$habitats
                if (sum(!(pvec%in%pops))>0) stop ("you have specified populations that can not occur in this landscape")
                Rland$individuals <- Rland$individuals[landscape.populations(Rland)%in%pvec,]
                if (!is.null(ns))
                {
                    if (length(ns)!=length(pvec))
                    {
                        ns=rep(ns[1],length(pvec))
                        message("created a ns vector out of a single popsize")
                    }
                    names(ns)=pvec
                    for (p in pvec)
                    {
                        inds <- Rland$individuals[landscape.populations(Rland)==p,]
                        Rland$individuals <-  Rland$individuals[landscape.populations(Rland)!=p,]
                        inddim <- dim(inds)[1]
                        if (!is.null(inddim))
                        {
                            if (inddim<ns[p])
                                Rland$individuals <- rbind(Rland$individuals,inds)
                            else
                                Rland$individuals <- rbind(Rland$individuals,inds[sample(1:inddim,ns[p],replace=F),])
                        }
                        
                    }
                    Rland$individuals <- Rland$individuals[order(Rland$individuals[,1],Rland$individuals[,4]),]
                }
                Rland
            }
    }
