# Crassostrea gigas population genomics

## The following are gigas population genomics analyses

### mtDNA 
data/NADH+NCR_fromNCBI.fas = two mitochondrial loci from whole mitochondrial genomes of 11 Crassostrea species at GenBank  
R/Blast_tree.R = Make NJ phylogeny  
output/Blast_tree.pdf = NJ phylogeny  

## The following are associated with the ABC analyses  
## Global grid (output/df.globe.Rda)  
lat = latitudinal midpoint of global 1º x 1º grid  
lon = longitudinal midpoint of global 1º x 1º grid  
coastline = Is there a coastline inside quadrat? (1 or 0)  
grid.for.model = Does our model include this quadrat? (1 or 0)  
reg = region (1_Asia, 2_wNA, 3_Eur)  
gridID = quadrat ID from global 1º x 1º grid (1 to 38750)  
plotted.Var1 and plotted.Var2 = Model grid position when plotted on a x-y figure    

## gigas source populations  (output/df.globe_source.csv)  
lat,lon,coastline,grid.for.model,reg,gridID,plotted.Var1,plotted.Var2 = same as df.globe    
sourceID = four putative source regions in Japan (hok, hon, sea, tok)  

## 1) C gigas SNPs - unpublished  
cgiga_meta.csv = metadata of grid    
cgiga_usat_Fst.Rda = pairwise Fst    
cgiga_usat_Fst3region.Rda = fst with Native vs wNA vs Eur    
cgiga_usat_FstNatNon.Rda = fst with Native vs Non-native    
cgiga_usat_FstOverall.Rda = overall Fst      
cgiga_usat_Hexp.Rda = Exp Heterozygosity    

## 2) G vermiculophylla microsatellites - Krueger-Hatfield 2017 E&E
gverm_meta.csv   
gverm_usat_Fst.Rda  
gverm_usat_FstOverall.Rda  
gverm_usat_FstNatNon.Rda  
gverm_usat_Fst3region.Rda  
gverm_usat_Hexp.Rda  

## 3) Undaria pinnitifida mtDNA  -Uwai 2006 Phycologia  
upinn_meta.csv  
upinn_mtdna_Fst.Rda  
upinn_mtdna_Fst3region.Rda  
upinn_mtdna_FstNatNon.Rda  
upinn_mtdna_FstOverall.Rda  
upinn_mtdna_nd.Rda = nucleotide diversity    

## 4) Sargassum muticum usat - LeCam et al 2020 Evolutionary Applications
smuti_meta.csv  
smuti_usat_Fst.Rda  
smuti_usat_Fst3region.Rda  
smuti_usat_FstNatNon.Rda  
smuti_usat_FstOverall.Rda  
smuti_usat_Hexp.Rda  

## 5) Batillaria attramentaria - Miura et al 2006 PNAS  
battr_meta.csv  
battr_mtdna_Fst.Rda  
battr_mtdna_Fst3region.Rda  
battr_mtdna_FstNatNon.Rda  
battr_mtdna_FstOverall.Rda  
battr_mtdna_nd.Rda  






