#' Get NDVI values between day1 and day2 at location (lng, lat)
#'
#' @param lng Longitude of the coordinate of interest.
#' @param lat Latitude of the coordinate of interest.
#' @param day1 Start day of the period.
#' @param day2 End day of the period.
#' @return  Vector with dates and NDVI values.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

get_NDVI_values<-function(lng, lat, day1, day2){
  sp <- SpatialPoints(cbind(lng, lat))
  pref<-"AH_NDVI_16d_250m_KE_"
  setwd(paste0(dir2,"/NDVI"))
  wh1<-max(which(datesNDVI<=day1))
  wh2<-min(which(datesNDVI>=day2))
  ans<-data.frame(date=c(), NDVI=c())
  for(whD in wh1:wh2){
    fln<-paste0(pref, datesNDVI[whD],".tif")
    vg <- raster(x=fln)
    rasValue <- raster::extract(vg, sp)
    rasValue<-rasValue/10000
    ans<-rbind(ans, data.frame(date=as.Date(datesNDVI[whD]),
                               NDVI=rasValue))
  }
  for(i in 2:nrow(ans)){
    cfs<-fit_line(ans$date[i-1], ans$date[i], ans$NDVI[i-1],
                  ans$NDVI[i])
    ds<-seq(as.Date(ans$date[i-1]),as.Date(ans$date[i]), by=1)
    for(j in 2:(length(ds)-1)){
      val<-cfs$intercept+cfs$slope*as.numeric(ds[j]-ds[1])
      ans<-rbind(ans, data.frame(date=as.Date(ds[j]),
                                 NDVI=val))
    }
  }
  ans<-ans[order(ans$date),]
  ans<-ans[which(ans$date>=day1 & ans$date<=day2),]
  return(ans)
}
