#' Temperature dependant duration of hopper development
#'
#' @param date Date of hopper hatching.
#' @param lng Longitude.
#' @param lat Latitude.
#' @return  Length of hopper development period in days
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
#'
hopper_dev_dur<-function(date, lng, lat){
  dev_day<-NA
  sp <- SpatialPoints(cbind(lng, lat))
  whD<-max(which(datesMET<=date))
  setwd(paste0(dir2,"/met_data"))
  fln<-paste0("Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
                                 collapse = ""),".nc")
  dat = nc_open(fln)
  tt  <- ncvar_get(dat, "time") 
  tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  xx<-ncvar_get(dat, "longitude")
  yy<-ncvar_get(dat, "latitude")
  varT <- ncvar_get(dat, "var__temperature__c_") 
  nc_close(dat)
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whDD<-which(as.Date(tms)>=date)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  ans<-data.frame(time=tms[whDD], tempr=varT[whx, why,whDD])
  dvlp_rt<-3*hopper_development(ans$tempr)/24 # 3-hour development rate
  cm_dvlp_rt<-cumsum(dvlp_rt)
  wh<-which(cm_dvlp_rt>=100)
  if(length(wh)>0){
    wh<-min(wh)
    dev_day<-as.Date(ans$time[wh])
  } else{
    if(whD<length(datesMET)){
      fln<-paste0("Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+1]),"-")),
                                     collapse = ""),".nc")
      dat = nc_open(fln)
      tt  <- ncvar_get(dat, "time") 
      tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
      xx<-ncvar_get(dat, "longitude")
      yy<-ncvar_get(dat, "latitude")
      varT <- ncvar_get(dat, "var__temperature__c_") 
      nc_close(dat)
      whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
      why<-which((yy-lat)^2==min((yy-lat)^2))[1]
      whDD<-which(as.Date(tms)>=date)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
      ans<-rbind(ans,data.frame(time=tms[whDD], tempr=varT[whx, why,whDD]))
      dvlp_rt<-3*hopper_development(ans$tempr)/24 # 3-hour development rate
      cm_dvlp_rt<-cumsum(dvlp_rt)
      wh<-which(cm_dvlp_rt>=100)
      if(length(wh)>0){
        wh<-min(wh)
        dev_day<-as.Date(ans$time[wh])
      } else {
        if(whD<length(datesMET)-1){
          fln<-paste0("Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+2]),"-")),
                                         collapse = ""),".nc")
          dat = nc_open(fln)
          tt  <- ncvar_get(dat, "time") 
          tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
          xx<-ncvar_get(dat, "longitude")
          yy<-ncvar_get(dat, "latitude")
          varT <- ncvar_get(dat, "var__temperature__c_") 
          nc_close(dat)
          whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
          why<-which((yy-lat)^2==min((yy-lat)^2))[1]
          whDD<-which(as.Date(tms)>=date)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
          ans<-rbind(ans,data.frame(time=tms[whDD], tempr=varT[whx, why,whDD]))
          dvlp_rt<-3*hopper_development(ans$tempr)/24 # 3-hour development rate
          cm_dvlp_rt<-cumsum(dvlp_rt)
          wh<-which(cm_dvlp_rt>=100)
          if(length(wh)>0){
            wh<-min(wh)
            dev_day<-as.Date(ans$time[wh])
          }
        }
      }
    }
  }
  return(dev_day)
}
