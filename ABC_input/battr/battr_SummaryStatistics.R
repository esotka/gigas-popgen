### Summary Statistics
library(strataG); library(ape)
rm(list=ls())
gt = readRDS("ABC_input/battr/battr_mtDNA_gtype")
seqs = read.FASTA("ABC_input/battr/Battr_180inds.fas")

# overall
overallHet = heterozygosity(gt, type = "expected") # haplotypic diversity

#### by pop ####
# 1) exp het (haplotypic diversity)
popHet = heterozygosity(gt, type = "expected",by.strata=T) 
# 2) overall PhiST
popOverall = overallTest(gt,nrep=0)$result["PHIst",1]
# 3) Nei pi = nucleotide diversity (within Pop)
nucD <- nucleotideDivergence(gt)
nucD_pop = nucD$within$mean; names(nucD_pop) = nucD$within$stratum
# 4) Nei nucleotide divergence (between Pop)
nucD_popPair = nucD$between$mean; names(nucD_popPair) = paste(nucD$between$strata.1,nucD$between$strata.2,sep="_")
# 5) Nei dA between pop
neiDA_popPair = nucD$between$dA; names(neiDA_popPair) = paste(nucD$between$strata.1,nucD$between$strata.2,sep="_")
# 6) pairwise PhiST
pws <- pairwiseSummary(pairwiseTest(gt,quietly=T,nrep=0))
phiST_popPair = pws$PHIst; names(phiST_popPair) = paste(pws$strata.1,pws$strata.2,sep="_")
# 7) Tajima D - on sequences
pop = factor(getStrata(gt))
tajD = c()
for(i in 1:length(levels(pop)))
{tajD = c(tajD,tajimasD(seqs[pop%in%levels(pop)[i]])$D)}
names(tajD) = levels(pop)

# by reg
meta = read.csv("ABC_input/battr/Miura2006haplotypes_meta_wSourceIDs.csv")
region = meta$sourceID[match(getStrata(gt),meta$Site)]
names(region) = names(getStrata(gt))
gt2 = gt
setStrata(gt2) = region
gt2 = gt2[!is.na(region)]
# 1) exp het (haplotypic diversity)
popHet_reg = heterozygosity(gt2, type = "expected",by.strata=T) 
# 2) overall PhiST
popOverall_reg = overallTest(gt2,nrep=0)$result["PHIst",1]
# 3) Nei pi = nucleotide diversity (within Pop)
nucD <- nucleotideDivergence(gt2)
nucD_reg = nucD$within$mean; names(nucD_reg) = nucD$within$stratum
# 4) Nei nucleotide divergence (between Pop)
nucD_regPair = nucD$between$mean; names(nucD_regPair) = paste(nucD$between$strata.1,nucD$between$strata.2,sep="_")
# 5) Nei dA between pop
neiDA_regPair = nucD$between$dA; names(neiDA_regPair) = paste(nucD$between$strata.1,nucD$between$strata.2,sep="_")
# 6) pairwise PhiST
pws <- pairwiseSummary(pairwiseTest(gt2,quietly=T,nrep=0))
phiST_regPair = pws$PHIst; names(phiST_regPair) = paste(pws$strata.1,pws$strata.2,sep="_")
# 7) Tajima D - on sequences
region = factor(region)
tajD_reg = c()
for(i in 1:length(levels(region)))
{tajD_reg = c(tajD_reg,tajimasD(seqs[region%in%levels(region)[i]])$D)}
names(tajD_reg) = levels(region)


print(c(
  popHet,
  popOverall,
  nucD_pop,
  neiDA_popPair,
  phiST_popPair,
  tajD,
  popHet_reg,
  popOverall_reg,
  nucD_reg,
  nucD_regPair,
  neiDA_regPair,
  phiST_regPair,
  tajD_reg))
  
  
  
