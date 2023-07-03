#' Get land cover type at location (lng, lat)
#'
#' @param lng Longitude of the coordinate of interest.
#' @param lat Latitude of the coordinate of interest.
#' @return Land cover type id.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

get_land_cover<-function(lng, lat){
  sp <- SpatialPoints(cbind(lng, lat))
  setwd(dir1)
  fl<-"PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
  lncv <- raster(x = paste0(fl))
  rasValue <- raster::extract(lncv, sp)
  rm(lncv)
  return(rasValue)
}