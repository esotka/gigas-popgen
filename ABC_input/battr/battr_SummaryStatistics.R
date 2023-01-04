### Summary Statistics
library(strataG)
rm(list=ls())
gt = readRDS("ABC_input/battr/battr_mtDNA_gtype")

# by pop
nucD <- nucleotideDivergence(gt)
print(nucD$within)
print(nucD$between)
pws <- pairwiseTest(gt)
pws.matrix <- pairwiseMatrix(pws,"PHIst")

# by reg
meta = read.csv("Miura2006haplotypes_meta_wSourceIDs.csv")
region = meta$sourceID[match(getStrata(gt),meta$Site)]
names(region) = names(getStrata(gt))
gt2 = gt
setStrata(gt2) = region
nucD <- nucleotideDivergence(gt2)
print(nucD$within)
print(nucD$between)
pws <- pairwiseTest(gt2)
pws.matrix <- pairwiseMatrix(pws,"PHIst")
print(pws.matrix)
