library(raster)
library(ggplot2)
library(ncdf4)
library(RNetCDF)
source('Utilities.R')

src<-read.csv("data/sourcelist_3860.csv", header = FALSE)
colnames(src)<-c("ID", "lng", "lat","Country")

ans<-data.frame(Year=c(), lng=c(), lat=c(), Length=c())
date<-as.Date("2020-10-01"); 
for(ii in 1:nrow(src)){
  lng<-src$lng[ii]; lat<-src$lat[ii]
  hoppdevday<-hopper_dev_dur(date, lng, lat, 24, 95)
  if(is.Date(hoppdevday)){
    ll<-as.numeric(hoppdevday-date)
    ans<-rbind(ans, data.frame(Year=2020, lng=lng, lat=lat, Length=ll))
  } else{
    ans<-rbind(ans, data.frame(Year=2020, lng=lng, lat=lat, Length=NA))
    
  }
  cat(c(ii, ""))
}

ggplot() +
  geom_tile(data = ans,
            aes(x = lng, y = lat, fill=Length)) +
  scale_fill_continuous(type = "viridis", direction=-1, na.value = 'gray') +
  labs(fill="DP (days)")+
  coord_fixed() +theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(ans[!is.na(ans$Length),], aes(x=Length))+
  geom_histogram(bins=20)+
  xlab("Developement period (days)")

mean(ans$Length[!is.na(ans$Length)])

