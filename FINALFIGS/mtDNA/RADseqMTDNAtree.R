rm(list=ls())
filenames <- list.files("FINALFIGS/mtDNA/fasta_fromAllan/")
library(ape)
system(paste("cat ",toString(paste("FINALFIGS/mtDNA/fasta_fromAllan/",filenames,sep="",collapse = " "))," > FINALFIGS/mtDNA/all.fasta",sep=""))

#dna <- read.FASTA("all.fasta")
dna <- read.dna("FINALFIGS/mtDNA/all.fasta","f")
tmp <- matrix(NA,nrow=717,4)
rownames(tmp) <- labels(dna)
colnames(tmp) <- c("A","C","G","T")
for (i in 1:717)
{tmp[i,] <- base.freq(dna[i,])}

# make visualization of base frequencies
matplot(tmp,type="l",col=1,xlab="ind",ylab="Base frequency")
legend(0,0.8,c("A","C","G","T"),lty=1:4,bty="n")

# remove outliers (any bp freq < 10%)
tmp2 <- tmp[complete.cases(tmp),]
tmp2 <- as.data.frame(tmp2)
outlier <- rowSums(tmp2<0.1)
tmp2 <- tmp2[outlier==0,]
matplot(tmp2,type="l",col=1,xlab="ind",ylab="Base frequency",main="no NAs and no freq < 10%")
legend(0,0.8,c("A","C","G","T"),lty=1:4,bty="n")
dna2 <- dna[labels(dna)%in%rownames(tmp2),]

### remove loci with few data
tmp <- c()
for (i in 1:dim(dna2)[2])
{tmp[i] <- table(as.character(dna2[,i])=="n" )["FALSE"]}
print(table(tmp)) ## loci with number of individuals genotyped

dna3 <- dna2[,tmp >= 615] ### 170 bp have at least 615 individuals

# remove individuals with any Ns
tmp <- c()
for (i in 1:dim(dna3)[1])
{tmp[i] <- table(as.character(dna3[i,])=="n" )["FALSE"]}
print(table(tmp)) ## 606 individuals have 170 bp

dna4 <- dna3[tmp==170 & !is.na(tmp),]

tmp <- matrix(NA,nrow=dim(dna4)[1],4)
rownames(tmp) <- labels(dna4)
colnames(tmp) <- c("A","C","G","T")
for (i in 1:dim(dna4)[1])
{tmp[i,] <- base.freq(dna4[i,])}

matplot(tmp,type="l",col=1,xlab="ind",ylab="Base frequency",main="606 ind; 170 bp")
legend(0,0.8,c("A","C","G","T"),lty=1:4,bty="n")

write.FASTA(dna4,"FINALFIGS/mtDNA/ind606.170bp.mtDNA.fas")


## find unique seqs and count
x <- readLines("FINALFIGS/mtDNA/ind606.170bp.mtDNA.fas")
d <- data.frame(names = x[c(TRUE, FALSE)], seq = x[c(FALSE, TRUE)])
out <- aggregate(names ~ seq, d, FUN = function(i)c(cnt = length(i), all = toString(i)))
write.table(out[,2],"FINALFIGS/mtDNA/ind606.170bp.mtDNA.uniq.stats.txt")
fasOUT <- data.frame(seq=c(paste("> Hap",1:length(out$seq),sep="_"),out$seq),order=c(seq(1,25,2),seq(2,26,2)))
fasOUT <- fasOUT[order(fasOUT$order),]
writeLines(fasOUT$seq,"FINALFIGS/mtDNA/uniq.170bp.fas")

### make NJ tree, combined with the concatenated NADH and NCR from genbank
library(muscle)
dna <- readDNAStringSet("FINALFIGS/mtDNA/NADH+NCR_concatenated_Final+UniqHaps.fas", "fasta")
print(dna)
aln <- muscle::muscle(dna)
aln.forApe <- as.DNAbin(aln)
## 1000 bootstrap reps
pdf("FINALFIGS/mtDNA/NADH+NCR_concatenated_Final+UniqHaps.pdf",height=5,width=4)
par(mar=c(0,0,0,0))
f <- function(x) nj(dist.dna(x))
tw <- f(aln.forApe) # NJ tree with K80 distance
set.seed(1)
## bootstrap with 100 replications:
(bp <- boot.phylo(tw, aln.forApe, f, quiet = TRUE,B = 1000))
bp[bp<900] = ""
plot(root(tw,"bilineata_MT985154.1"),cex=0.5)
drawSupportOnEdges(100*round(as.numeric(bp)/1000,2),cex=0.5,bg="white",frame="none",adj=c(0.5,-.75),col="red")
segments(0.05,5,0.07,5)
text(0.06,6,"0.02")
dev.off()


