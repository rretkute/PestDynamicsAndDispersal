#' Hopper development conditions at location (lng, lat)
#' .
#'
#' @param date Date hoppers hatched from eggs.
#' @param lng Longitude.
#' @param lat Latitude.
#' @param thresh Minimal NDVI required for hopper development.
#' @param peak If account for appearance of fresh vegetation (default is TRUE).
#' @param lambda Parameter for smoothing NDVI.
#' @param m Number of change points in NDVI profile.
#' @return Development date (NA if unseccesfull).
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

hopper_dev<-function(date, lng, lat, thresh, peak=TRUE,
                     lambda=1, m=2){
  ans<-NA
  dev.period<-hopper_dev_dur(date, lng, lat)
  if(!is.na(dev.period)){
    day1<-date
    day2<-day1+dev.period
    ndvi<-get_NDVI_values(lng, lat, day1, day2)
    if(all(ndvi$NDVI>=thresh)){
      ans<-day2
    }  else {
    if(peak==TRUE){
      ndvi<-get_NDVI_values(lng, lat, min(datesMET), max(datesNDVI))
      yy<-smooth_Whittaker(ndvi$NDVI, lambda=1)
      p<-find_peaks(yy, m=2)
      whD<-which(ndvi$date[p]>=day1 & ndvi$date[p]<=day2)
      if(length(whD)>0) ans<-day2
      }
    }
  }
  return(ans)
}
