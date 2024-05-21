##
### plot ship tracks
###

library(maps)
library(SDMTools)
library(sf)
library(dplyr)
library(parallel)
cores=10
direction=T
portlength = 2  # number of days between entries to signify the start of a leg
maxdeltalong = 40#change in longitude to cause an NA in plotlong
maxdeltalat = 40#change in lat to cause an NA in plotlat
maxVoyageDays = c(30,60)
years = c(1950:2010)
japanbox=list(long=c(129.44,146.12),lat=c(30.6,46.144))
origin_files = c("hok.kml", "tokyo.kml", "sea.kml", "miyagi.kml","kag.kml")
#origin_files = c("tokyo.kml")
dest_files = c("NZ.kml","EU.kml","pnw.kml","SoCal.kml")

###shiptracks = readRDS("shiptracks.RDS")
shiptracks = do.call(rbind,lapply(1:7,function(i) readRDS(paste0("shiptracks_",i,".RDS"))))
tracks = shiptracks %>% filter(year %in% years) %>% arrange(year,ID,time) %>% filter(!grepl("\\*",ID)) %>% mutate(long360=ifelse(long<0,360+long,long))

IDyear = tracks %>% filter(long360<=japanbox$long[2],long360>=japanbox$long[1])%>%
    filter(lat<=japanbox$lat[2],lat>=japanbox$lat[1]) %>% select(ID,year) %>% distinct()

tracks = tracks%>%right_join(IDyear)

tracks$deltalong = c(NA,diff(tracks$long360))
tracks$laglong = c(NA,tracks$long360[-nrow(tracks)])
tracks$plotlong = tracks$long360
tracks$plotlong[abs(tracks$deltalong)>=maxdeltalong] = NA

tracks$deltalat = c(NA,diff(tracks$lat))
tracks$laglat = c(NA,tracks$lat[-nrow(tracks)])
tracks$plotlat = tracks$lat
tracks$plotlat[abs(tracks$deltalat)>=maxdeltalat] = NA


tracks.orig=tracks


if (FALSE) #plot tracks, set to FALSE to avoid the overhead
{
    pdf("plot_japan_tracks.pdf")
    for (y in c(1958,1968,1978)) #unique(tracks$year))
    {
        print(y)
        map('world2')

        ydf=tracks[tracks$year==y,]
        lapply(unique(ydf$ID), function(id)
        {
            pldf=ydf[ydf$ID==id,]
            pldf=pldf[order(pldf$time),]
            points(y=pldf$plotlat,x=pldf$plotlong,type="l",lwd=0.5,col="red")
        })
    }
    dev.off()
}

trackll=as.matrix(tracks[,c("long","lat")])

origins = do.call(cbind,mclapply(origin_files,mc.cores=cores,function(f)
{
    tmp=as.matrix(st_coordinates(st_read(f,quiet=T))[,1:2])
    rv = data.frame(v1=pnt.in.poly(trackll,tmp)[,3])
    names(rv)=gsub(".kml","",f)
    rv
}))


dests = do.call(cbind,mclapply(dest_files,mc.cores=cores,function(f)
{
    print(f)
    tmp=as.matrix(st_coordinates(st_read(f,quiet=T))[,1:2])
    rv = data.frame(v1=pnt.in.poly(trackll,tmp)[,3])
    names(rv)=gsub(".kml","",f)
    rv
}))

names(origins)=paste0("org.",names(origins))
names(dests)=paste0("dest.",names(dests))


##this set of loops is intended to identify id/year combos that have both a particular
## origin and a particualr destination

id_yrs_orig = NULL
for (i in names(origins))
    {
        tdf  = unique(tracks[origins[,i]>0,c("ID","year")])
        tdf$origin=i
        id_yrs_orig=rbind(id_yrs_orig,tdf)
    }

id_yrs_dest = NULL
for (i in names(dests))
    {
        tdf  = unique(tracks[dests[,i]>0,c("ID","year")])
        tdf$dest=i
        id_yrs_dest=rbind(id_yrs_dest,tdf)
    }

orig_dest = merge(id_yrs_orig,id_yrs_dest,all=F) #same ID went to origin and dest
orig_dest$comb = factor(paste0(orig_dest$origin,"_",orig_dest$dest))

if (TRUE) #plot tracks, set to FALSE to avoid the overhead
{
    pdf("plot_org_dest_tracks.pdf")
    for (y in 1950:1965) #unique(tracks$year))
    {
        print(y)
        map('world2')
        ydf=tracks[tracks$year==y,]
        ydf=ydf[ydf$ID%in%orig_dest$ID,]
        if (nrow(ydf)>0)
            lapply(unique(ydf$ID), function(id)
            {
                pldf=ydf[ydf$ID==id,]
                pldf=pldf[order(pldf$time),]
                for (pair in unique(unclass(orig_dest$comb)[orig_dest$ID==id]))
                {
                    print(pair)
                    print(class(pair))
                    points(y=pldf$plotlat,x=pldf$plotlong,type="l",lwd=0.5,col=1+pair)
                }
                    

            })
    }
    dev.off()
}



tracks=cbind(tracks,origins,dests)
tracks=tracks%>%filter(ID%in%orig_dest$ID) %>% arrange(ID,time)

id = unique(tracks$ID)[1]
####track times do it for each ID
ttimes=do.call(rbind,
               mclapply(unique(tracks$ID),mc.cores=cores,
                      function(id)
                      {
                          strk = tracks%>% filter(ID == id)
                          if (length(grep("org",names(strk)))>1)
                              strk$orgs=rowSums(strk[,names(strk)[grep("org",names(strk))]])>0 else strk$orgs=strk[,names(strk)[grep("org",names(strk))]]>0
                          if (length(grep("dest",names(strk)))>1)
                              strk$dests=rowSums(strk[,names(strk)[grep("dest",names(strk))]])>0 else strk$dests=strk[,names(strk)[grep("dest",names(strk))]]>0
                          
                          strk = strk %>% filter(orgs|dests)
                          tt = rbind(strk %>% group_by(orgs) %>% arrange(time) %>% filter(row_number()==n()),    #last obs in source reg
                                     strk %>% group_by(dests) %>% arrange(time) %>% filter(row_number()==1)) %>% #first obs in dest region
                              arrange(time)
 #                         tt$tdiff = c(diff(tt$time),0)

                          rdf=NULL
                          
 #                         tt = tt[tt$tdiff<=maxVoyageDays,]

                          if (nrow(tt)>1)
                              {
                                  tt$order = with(tt,ifelse(orgs & (!dests),1,ifelse((!orgs)&dests,2,NA)))
                                  for (rw in 1:nrow(tt))
                                  {
                                      keep=TRUE
                                      if (rw==nrow(tt)) #last row in df
                                      {
                                          if (!tt$dests[rw]) keep=FALSE 
                                      } else #all the rest of the rows
                                      {
                                          if (tt$orgs[rw])
                                              if (tt$dests[rw+1])
                                              {
                                                  ds = names(tt)[grep("dest\\.",names(tt))][which(tt[rw+1,grep("dest\\.",names(tt))]>0)]
                                                  os = names(tt)[grep("org\\.",names(tt))][which(tt[rw,grep("org\\.",names(tt))]>0)] 
                                                  rdf=rbind(rdf,data.frame(ID=id,otime=tt$time[rw],dtime=tt$time[rw+1],
                                                                           olat=tt$lat[rw],olong=tt$long[rw],
                                                                           dlat=tt$lat[rw+1],dlong=tt$long[rw+1],
                                                                           o=os[1],d=ds[1]
                                                                           )
                                                            )
                                              }
                                      }
                                  }
                              }
                          if (!is.null(rdf))
                          {
                              rdf$timediff=rdf$dtime-rdf$otime
                              rdf[rdf$timediff>0,]
                          } else rdf
                      }))

ttimes=ttimes[ttimes$timediff<=max(maxVoyageDays),]

screened = do.call(rbind,
                   mclapply(1:nrow(ttimes),mc.cores=cores,
                            function(l)
                            {
                                tracks %>% filter(ID==ttimes$ID[l],time>=ttimes$otime[l], time<=ttimes$dtime[l] )
                            }))

###lets look at some tracks
if (FALSE) #commented out for now
{
    for(y in years)
    {
        pdf(paste0("tracksfor_",y,".pdf"))
        strk = screened[screened$year==y,]
        map("world")
        for(id in unique(strk$ID))
        {
            points(lat~plotlong,data=strk[strk$ID==id,],type="l",lwd=0.7,col=which(unique(strk$ID)==id))
        }
        dev.off()
    }
}


tracks.bak=tracks

routes = expand.grid(o=names(origins),d=names(dests))
routes$route =  apply(routes,1,function(x){paste0(x[1],"->",x[2])})

ttimes = left_join(ttimes,routes)

for (days.larvae.survive in maxVoyageDays)
{
print(days.larvae.survive)
tbl = with(ttimes[ttimes$timediff<=days.larvae.survive,],table(o,d))

colnames(tbl) = gsub("dest.","",colnames(tbl))
rownames(tbl) = gsub("org.","",rownames(tbl))
        
save(file=paste0("days_",days.larvae.survive,".rda"),tbl,days.larvae.survive,tracks,screened, ttimes)
dump("tbl",file=paste0("shippingTbl",days.larvae.survive,".R"))
print("about to plot")
pdf(paste0("ship_log_trips_len_",days.larvae.survive,".pdf"))

barplot(tbl,beside=T,col=heat.colors(5),main=paste("numbers of trips",days.larvae.survive,"days or less"),xlab="Destination")
legend(y=0.95*max(tbl),x=3,legend=rownames(tbl),fill=heat.colors(5)[1:5],title="Sources")

dev.off()
}

