### to run on hpc.cofc.edu
# module load apptainer
# singularity run /home/stranda/bin/bioos.simg R

library(strataG)
gt = readRDS("battr_mtDNA_gtype")
meta = read.csv("battr_meta.csv")
reg = meta$source[match(getStrata(gt),meta$pop)
#Cali       hon       kag nonSource       PNW       sea       tok 
#  20        60        10        10        20        40        20 
intro = gt[reg%in%c("Cali","PNW"),]


