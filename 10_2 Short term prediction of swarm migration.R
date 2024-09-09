library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(RColorBrewer)
library(ncdf4)
library(RNetCDF)
library(geosphere)
source('Utilities.R')

src<-read.csv("data/sourcelist_3860.csv", header = FALSE)
colnames(src)<-c("ID", "lng", "lat","Country")

pth<-read.csv("data/forecast_path.csv", header = TRUE)
lnd.pnts<-pth[c(1, 17),]
obs<-pth[6,]

fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=pth, size=1,
            aes(x=lng, y=lat), col='black')+
  geom_point(data=obs, aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  geom_point(data=lnd.pnts[1,], aes(x=lng, y=lat), fill='lightblue',
             pch=21, size=4)+
  geom_point(data=lnd.pnts[2,], aes(x=lng, y=lat), fill='darkblue',
             pch=21, size=4)+
  xlab("")+ ylab("") +
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(45, 50),  ylim = c(6, 11), expand = FALSE)

n.sim<-100000
ans<-data.frame(lng=c(), lat=c(), dist1=c(), dist2=c())
set.seed(1)
while(nrow(ans)<n.sim){
  date<-as.Date("2020-11-19")
  lng<-runif(1, min=obs$lng-5, max=obs$lng+5)
  lat<-runif(1, min=obs$lat-5, max=obs$lat+5)
  yy<-distm(c(lng, lat), obs[,c('lng', 'lat')], fun = distHaversine)/1000
  if(yy<=500){
    day.trj<-get_wind(date, lng, lat, c())
    xx<-min(distm(day.trj[, c('lng', 'lat')], obs[,c('lng', 'lat')], fun = distHaversine)/1000)
    ans<-rbind(ans, data.frame(lng=lng, lat=lat, dist1=yy[1], dist2=xx))
    cat(c(nrow(ans), ""))
  }
}

dst<-data.frame(Distance=c(), Fraction=c())
for(dd in seq(10, 500)){
  wh<-which(ans$dist1<=dd & ans$dist2<=5)
  dst<-rbind(dst, data.frame(Distance=dd, Fraction=length(wh)/nrow(ans)))
}

ggplot(dst, aes(x=Distance, y=Fraction))+
  geom_path() +
  scale_y_continuous(trans='log10')+
  xlab("Radius (km)")+
  ylab("Fraction within 5km")

n.sim<-10000
TRJ<-data.frame(date =c(), lng=c(), lat=c(), flt.tm=c(), fl.ds=c(),
                 id=c(), day=c(), obs=c())
set.seed(1)
id<-0
while(id<=n.sim){
  date<-as.Date("2020-11-19")
  lng<-runif(1, min=obs$lng-2, max=obs$lng+2)
  lat<-runif(1, min=obs$lat-2, max=obs$lat+2)
  yy<-distm(c(lng, lat), obs[,c('lng', 'lat')], fun = distHaversine)/1000
  if(yy[1]<=250){
  day.trj<-get_wind(date, lng, lat, c())
    id<-id+1
    fl.ds<-1
    tm<-as.numeric(day.trj$date-min(day.trj$date))/(60*60)
    day.trj$flt.tm<-tm
    day.trj$id<-id
    day.trj$day<-date
    day.trj$fl.ds<-fl.ds
    day.trj$obs<-0
    xx<-distm(day.trj[,c('lng', 'lat')], obs[,c('lng', 'lat')], fun = distHaversine)/1000
    if(min(xx)<=5) day.trj$obs<-1
    TRJ<-rbind(TRJ, day.trj)
    cat(c(id, ""))
  }
}

fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=TRJ, size=0.25,
            aes(x=lng, y=lat, group=id, col=flt.tm), col="gray", alpha=0.2) +
  geom_path(data=pth, size=1,
            aes(x=lng, y=lat, group=id), col='black')+
  geom_point(data=obs, aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  xlab("")+ ylab("") +
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(45, 50),  ylim = c(6, 11), expand = FALSE)

fig+
  annotation_scale(location = "br", width_hint = 0.5)+
  geom_path(data=TRJ[TRJ$obs==1,], size=0.25,
            aes(x=lng, y=lat, group=id, col=flt.tm)) +
  scale_color_viridis(option = "D", direction=-1)+
  geom_path(data=pth, size=1,
            aes(x=lng, y=lat, group=id), col='black')+
  geom_point(data=obs, aes(x=lng, y=lat), pch=21,
             fill="red", size=4)+
  xlab("")+ ylab("") +
  guides(col="none", fill="none")  +
  coord_sf(xlim = c(45, 50),  ylim = c(6, 11), expand = FALSE)

