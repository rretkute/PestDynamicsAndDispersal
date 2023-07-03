#' Sample individual wind trajectory at day date and location (lng, lat)
#' .
#'
#' @param date Date hoppers hatched from eggs.
#' @param lng Longitude.
#' @param lat Latitude.
#' @param rn.sd Set random seed.
#' @return Development rate (per day).
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

get_wind<-function(date, lng, lat, rn.sd){
  setwd(paste0(dir2,"/wind_trajectories/"))
  tmp<-sqrt((src$lng-lng)^2 +(src$lat-lat)^2)
  wh<-which(tmp==min(tmp))
  d.lng<-src$lng[wh]-lng
  d.lat<-src$lat[wh]-lat
  if(wh<10) s<-paste0("0000",wh)
  if(wh>=10 & wh<100) s<-paste0("000",wh)
  if(wh>=100 & wh<1000) s<-paste0("00",wh)
  if(wh>=1000 & wh<10000) s<-paste0("0",wh)
  dt<-format(date,"%Y%m%d")
  flnm <- paste0("locust_traj_",dt,"_nc/","Data_Traj_Swarm", s, "_C1_",dt,".nc")
  dat = nc_open(flnm)
  tt <- ncvar_get(dat, "time")            # Extract time scale
  tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  yy<-ncvar_get(dat, "latitude")
  xx<-ncvar_get(dat, "longitude")
  nc_close(dat)
  if( length(rn.sd)==1) set.seed(rn.sd)
  j<-sample(1:1000,1)
  wn.trj<-data.frame(date=tms, lng=xx[j,]-d.lng, lat=yy[j,]-d.lat)
  wn.trj<-wn.trj[order(wn.trj$date),]
  return(wn.trj)
}