# gigas-popgen

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
