library(raster)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
source('Utilities.R')

grd.1km<-read.csv('data/grd_1km.csv', header = TRUE)

date<-as.Date("2021-01-15")

st.dur<-data.frame(ID=c(), date=c(), lng=c(), lat=c(), Duration=c())
for(jj in 1:nrow(grd.1km)){
  lng<-grd.1km$lng[jj] 
  lat<-grd.1km$lat[jj] 
  lct<-get_land_cover(lng, lat)
  ndvi<-get_NDVI_values(lng, lat, date-150, date)
  if(!is.na(ndvi$NDVI[1])){
    if(sum(ndvi$NDVI)>0 & min(ndvi$NDVI)>=0){
      yy<-smooth_Whittaker(ndvi$NDVI, lambda=1) 
      ndvi$NDVIsmooth<-yy
      a<-swrm_stay_type(lct, ndvi, 6) 
      if(!is.na(a)){
        st.dur<-rbind(st.dur,
                      data.frame(ID=jj, date=date, lng=lng, lat=lat, Duration=a))
      }
    }
  }
  cat(c(jj,""))
}

st.dur$Duration<-factor(st.dur$Duration, 
                        levels=c("Short stay", "Medium stay", "Long stay"))
fig+
  geom_tile(data=st.dur, aes(x=lng, y=lat, fill=Duration)) +
  xlab("")+ ylab("") +
  scale_fill_manual(name = "Duration",
                    values=c("#7570b3", "#1b9e77", "#d95f02", "gray"))+
  theme(legend.position="none") 


