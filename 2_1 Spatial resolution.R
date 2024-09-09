library(raster)
library(ggplot)

coord.lim<-c(30, 55, -5, 20)

grd.init <- raster(xmn= coord.lim[1], ymn= coord.lim[3], 
            xmx = coord.lim[2], ymx = coord.lim[4], 
            resolution = 0.008,
            crs = '+proj=utm +zone=38 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ')
grd.init <- setValues(grd.init, 1)
grd.init

cntr<-c("ET", "KE", "SO", "ER", "DJ")

for(ii in 1:length(cntr)){
  bnrs <- raster::getData('GADM', country= cntr[ii], level=0)
  rc <- crop(grd.init, extent(bnrs))
  rc <- mask(rc, bnrs)
  plot(rc)
  y = rasterToPoints(rc)
  y<-as.data.frame(y)
  colnames(y)<-c("lng", "lat", "Country")
  y$Country<-cntr[ii]
  if(ii==1){
    grd.1km<-y
  } else {
    grd.1km<-rbind(grd.1km, y)
  }
}
nrow(grd.1km)
head(grd.1km)

ggplot(grd.1km[sample(1:nrow(grd.1km), 2*10^5),], aes(x=lng, y=lat))+
  geom_point(aes(col=Country), pch=15, size=0.01) +
  theme(legend.position="bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))+
  coord_fixed()

ggplot(grd.1km, aes(x=lng, y=lat))+
  geom_point(aes(col=Country), pch=15, size=0.1)+
  xlim(41.75, 42.25) + ylim(3.75, 4.25) +
  theme(legend.position="none") +
  coord_fixed()

write.csv(grd.1km, file='data/grd_1km.csv', row.names = FALSE)

