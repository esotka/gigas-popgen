library(strataG)
ingen = read.csv("Zostera_japonica_SSR_rename_loci.csv")
print(names(ingen))
gt = df2gtypes(ingen,ploidy=2)

meta=read.csv("zostera_meta.csv")
