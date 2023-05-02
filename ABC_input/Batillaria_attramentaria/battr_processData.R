#### Download haplotypes and make a fasta file
## Note: frequency table graciously sent by Dr. Miura (June 2022)
## lat / lon for Japan sent by Dr. Miura. 
##  lat / lon for wNA generated using google maps
rm(list=ls())
library(ape)
library(readxl)
library(reshape)
#library(muscle) # BiocManager
metaHap <- read_xlsx("Miura2006haplotypes_haps.xlsx") 
loc1 <- read.GenBank(metaHap$Haplotype)
write.dna(loc1,"battr_mtdna.unique.fas",format = "fasta")
dist.aln <- dist.dna(loc1) # same length
plot(nj(dist.aln))
# prepare for strataG
# if you want to stratify the resulting gtypes object from sequence2gtypes(), one sequence for each individual should be provided, rather than just a set of unique haplotypes.
metaHap <- data.frame(metaHap)
md <- melt(metaHap,id="Haplotype")
md <- md[complete.cases(md),]

meta <- read_xlsx("Miura2006haplotypes_meta.xlsx")

popID <- rep(as.character(md$variable),md$value)
hapInd <- rep(md$Haplotype,md$value)
lonInd <- as.numeric(meta$longitude[match(popID,meta$Site)])
latInd <- as.numeric(meta$latitude[match(popID,meta$Site)])
indID <- paste("ind",1:length(hapInd),sep="_")
country <- ifelse(lonInd>0,"1_Asia","2_wNA")
metaInd <- data.frame(indID,popID,country,lonInd,latInd,hapInd)

indDNA <- loc1[match(metaInd$hapInd,names(loc1))]
names(indDNA) <- indID
write.dna(indDNA,"Battr_180inds.fas",format="fasta")
write.csv(metaInd,"Battr_180inds_meta.csv",quote=F,row.names=F)

##save as strataG::gtypes#
library(strataG)
haps <- read.dna("Battr_180inds.fas",format="fasta")
meta <- read.csv("Battr_180inds_meta.csv")
strata <- meta$popID
names(strata) <- meta$indID
loc1.g <- sequence2gtypes(haps,strata)
saveRDS(loc1.g,"battr_mtDNA_gtype")


