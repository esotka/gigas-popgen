library(strataG)
ingen = read.csv("Zostera_japonica_SSR_rename_loci.csv")
print(names(ingen))

gto = df2gtypes(ingen,ploidy=2)

###I went ahead and inputed assuming hwe in each pop and using only the alleles found in each pop
miss=data.frame(which(is.na(ingen),arr.ind = T))
miss$pop=ingen$pop[miss$row]
for (p in unique(ingen$pop))
    for (l in letters[1:25])
    {
        lcols= grep(l,names(ingen))
        sset = unname(unlist(ingen[ingen$pop==p,lcols]))
        sset=sset[!is.na(sset)]
        lmiss= miss[miss$pop==p,]
        if (nrow(lmiss)>0) #there are missing data in this pop
            if (any(lcols %in% lmiss$col)) #and there are missing data for this loc in this pop
                for (i in 1:nrow(lmiss))
                {
                    ingen[lmiss$row[i],lmiss$col[i]] = sample(sset,1)
                }
    }
gt = df2gtypes(ingen,ploidy=2)
meta=read.csv("zostera_meta.csv")
if (all(meta$pop %in% getStrata(gt))&all(getStrata(gt)%in%meta$pop))
    saveRDS(file="Zostera_japonica_SSR_gtype",gt)
