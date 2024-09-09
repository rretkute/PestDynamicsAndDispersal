library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(raster)
library(ggspatial)
library(sf)
sf_use_s2(FALSE)
source('utilities.R')

# Large file ~1.7Gb. It should be downloaded from the provider, see Supplement Information
flnm<-"PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
dir<-"Downloads/"
fl<-paste0(dir, flnm)
lncv <- raster(x = paste0(fl))
co_extent <- extent(coord.lim)
co_extent <- as(co_extent, "SpatialPolygons")
sp::proj4string(co_extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
proj_co_extent <- spTransform(co_extent, crs(data))
lncv <- crop(lncv, proj_co_extent)

grd.1km<-read.csv('data/grd_1km.csv', header = TRUE)
sp <- SpatialPoints(cbind(grd.1km$lng, grd.1km$lat))
rasValue <- raster::extract(lncv, sp)

y<-data.frame(lng=grd.1km$lng, lat=grd.1km$lat, Type=rasValue)
ys<-sort(unique(y$Type))

clsf<-read.csv("data/discrete_classification_Class_Table.csv",  stringsAsFactors=FALSE)
colnames(clsf)<-c("Value","Color","Description","Full_description")

cols <- setNames(clsf$Color, clsf$Description)

y$LandCoverType<-"Unknown"
for(ii in 1:nrow(clsf)){
  y$LandCoverType[y$Type==clsf$Value[ii]]<-clsf$Description[ii]
}

fig+
  geom_point(data = y, 
             aes(x = lng, y = lat,  col=LandCoverType), shape = 15, size=0.0001) +
  scale_color_manual(values = cols) +
  guides(color = guide_legend(override.aes = list(size = 3)))
 theme(legend.position="none")

# write.csv(y[, 1:3], file='data/Land_cover_1km.csv', row.names = FALSE)
