### pull in vcf files.
# 1) filter out bad individuals (<5% of SNPs are NAs)
# 2) keep loci that are called at all individuals

library(hierfstat)
library(vcfR)

if (!exists("filelist"))
{
    print("using built-in filelist")
    filelist <- c("all.bamlist_NC_047559.1.vcf.gz",
                  "all.bamlist_NC_047560.1.vcf.gz",
                  "all.bamlist_NC_047561.1.vcf.gz",
                  "all.bamlist_NC_047562.1.vcf.gz",
                  "all.bamlist_NC_047563.1.vcf.gz",
                  "all.bamlist_NC_047564.1.vcf.gz",
                  "all.bamlist_NC_047565.1.vcf.gz",
                  "all.bamlist_NC_047566.1.vcf.gz",
                  "all.bamlist_NC_047567.1.vcf.gz",
                  "all.bamlist_NC_047568.1.vcf.gz")
}

out <- c()
for (i in 1:length(filelist)) #
{
dat <- read.vcfR(filelist[i])
gt2 <- extract.gt(dat,return.alleles=F)
gt2[gt2=="0/0"] <- 0
gt2[gt2=="0/1"] <- 1
gt2[gt2=="1/1"] <- 2
out <- rbind(out,gt2)
}

#loci <- rownames(out)
inds <- colnames(out)#readLines("data/ind393_clean")
#inds2 <- substr(inds,52,nchar(inds))
#inds3 <- gsub("_bwa.np.bam","",inds2)
#colnames(out) <- inds3
inds <- gsub("_bwa.np.bam","",gsub("/home/stranda/projects/oyster/sotka/roslin/bwa/out/","",inds))
geno <- t(matrix(as.numeric(out),nrow=nrow(out)))
geno <- data.frame(geno)
rownames(geno) <- inds
colnames(geno) <- rownames(out)

### remove ind with >5% NAs
NAs.ind <- rowSums(is.na(geno)) # how many missing
geno2 <- geno[NAs.ind<=ncol(geno)*.05,]
pop <- substr(rownames(geno2),1,3)

### Keep loci that have calls
NAs.loc <- colSums(is.na(geno2))
geno2 <- geno2[,NAs.loc==0]

write.table(data.frame(id=rownames(geno2),geno2),
            file="SNPs_noZeros.txt",quote=F,sep="\t",row.names = F)




