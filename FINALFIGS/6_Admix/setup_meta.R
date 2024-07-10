###
### R script to setup metadata for pca and admixture analysis
### Uses a bamfile and metadata to classify points
###


bf = "all.bamlist"
metaf = "gigas_bam_qualities.csv"
meta2 = "../gigas_meta_41pop_env_FINAL.csv"
##read data and merge bams with metadata


bams=readLines(bf)

meta=read.csv(metaf)
meta2=read.csv(meta2)

meta=meta[,c("pop","Region","id","Latitude","Longitude")]
meta$Region2=meta2$Region2[match(meta$pop,meta2$pop)]
bams <- gsub("/home/stranda/projects/oyster/sotka/roslin/bwa/out/","",bams)
bams <- gsub("/home/stranda/projects/oyster/LiSamples/roslin/bwa/out/","",bams)
bams <- gsub("_roslin","",bams)
bams <- gsub("_bwa.np.bam","",bams)

iddf = data.frame(id=bams,bforder=1:length(bams))

print(names(iddf))

iddf <- merge(iddf,meta,all.x=T)
iddf <- iddf[order(iddf$bforder),]

print("assembled dataframe")

iddf$Region = factor(as.character(iddf$Region))
iddf$Region = factor(iddf$Region,
                     levels=levels(iddf$Region)[c(8,1,2,5,3,9,13,12,4,11,10,6,7)])
iddf$Region2 = factor(as.character(iddf$Region2))
iddf$Region2 = factor(iddf$Region2,
                     levels=levels(iddf$Region2)[c(8,4,10,13,5,1,11,9,12,7,2,3,6)])


iddf$Latitude[is.na(iddf$Latitude)]=0
iddf$Longitude[is.na(iddf$Longitude)]=0
print(unique(iddf$Region))

#now there's a dataframe called iddf that is pretty complete metadata
