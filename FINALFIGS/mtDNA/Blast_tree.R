rm(list=ls())
library(ape,quietly=T)
library(muscle,quietly = T)
pdf("output/Blast_tree.pdf",height=12,width=7)
par(mfrow=c(2,1),mar=c(1,1,1,1))
# NADH tree
dna <- readDNAStringSet("data/NADH+NCR_fromNCBI.fas", "fasta")
id <- strsplit(names(dna),"_")
locus <- unlist(lapply(id,"[[",3))
spp <- unlist(lapply(id,"[[",2))
dna2 <- dna[!spp=="plicatula"]
id2 <- strsplit(names(dna2),"_")
locus2 <- unlist(lapply(id2,"[[",3))
spp2 <- unlist(lapply(id2,"[[",2))
aln <- muscle::muscle(dna2[locus2=="NADH"])
aln.forApe <- as.DNAbin(aln)
dist.aln <- dist.dna(aln.forApe)
#print(dist.aln)
plot(nj(dist.aln),cex=.5,main="NADH",show.tip.label=F)
tiplabels(spp2[locus2=="NADH"],cex=.5,frame="none")

# non-coding region
aln <- muscle::muscle(dna2[locus2=="nonCodingRegion"])
aln.forApe <- as.DNAbin(aln)
dist.aln <- dist.dna(aln.forApe)
#print(dist.aln)
plot(nj(dist.aln),cex=.5,main="Non-coding region",show.tip.label=F)
tiplabels(spp2[locus2=="nonCodingRegion"],cex=.5,frame="none")

dev.off()

