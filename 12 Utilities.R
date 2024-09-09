
# Limits for maps etc
coord.lim<-c(30, 55, -5, 20)
co_extent <- extent(coord.lim)
co_extent <- as(co_extent, "SpatialPolygons")
sp::proj4string(co_extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


world <- ne_countries(scale = "medium", returnclass = "sf")
# Remove a border between Somalia ans Somaliland
sml <- world %>% 
  dplyr::select(sov_a3) %>% 
  filter(sov_a3 %in% c("SOL", "SOM")) %>% 
  mutate(sov_a3 = "SML") %>%  
  group_by(sov_a3) %>% 
  summarise()

world <- world %>% 
  dplyr::select(sov_a3) %>% 
  filter(!sov_a3 %in% c("SOL", "SOM")) %>% 
  bind_rows(sml)

fig<-ggplot(data = world) +
  geom_sf(fill="white") +
  annotation_scale(location = "br", width_hint = 0.5) +
  coord_sf(xlim = c(coord.lim[1], coord.lim[2]), 
           ylim = c(coord.lim[3], coord.lim[4]), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

datesMET<-c("2020-08-01", "2020-09-01", "2020-10-01", "2020-11-01", "2020-12-01")

pref<-"AH_NDVI_16d_1km_"
lf<-list.files("data/NDVI/")
dates<-c()
for(i in 1:length(lf)){
  x<-strsplit(lf[i],'\\.')
  x<-x[[1]][1]
  x<-strsplit(x,pref)
  x<-x[[1]][2]
  dates<-c(dates, x)
}
datesNDVI<-sort(dates)

suit_breeding<-function(date, lng, lat){
  sp <- SpatialPoints(cbind(lng, lat))
  doy <- as.numeric(strftime(date, format = "%j"))
  # Breeding suitability map
  fln<-"data/breeding_suitability_map_1km.tif"
  sbr <- raster(x=fln)
  sbr <- crop(sbr, co_extent)
  # plot(sbr)
  sbr<-raster::extract(sbr, sp, method='bilinear')/100
  if(is.na(sbr)) sbr<-0
  # 24 hours prior
  whD<-max(which(datesMET<=date-1))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
                                 collapse = ""),".nc")
  dat = nc_open(fln)
  tt  <- ncvar_get(dat, "time") 
  tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  xx<-ncvar_get(dat, "longitude")
  yy<-ncvar_get(dat, "latitude")
  varP <- ncvar_get(dat, "var__precipitation_rate__mm_hr_") 
  varSM <- ncvar_get(dat, "var__soil_moisture")
  nc_close(dat)
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whD<-which(as.Date(tms)==date-1)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  prec<-max(varP[whx, why, whD])
  whD<-which(as.Date(tms)==date-1)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  if(length(whD)==0) whD<-1
  soilm<-max(varSM[whx, why, whD])
  # 48 hours prior
  whD<-max(which(datesMET<=date-2))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
                                      collapse = ""),".nc")
  dat = nc_open(fln)
  tt  <- ncvar_get(dat, "time") 
  tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  xx<-ncvar_get(dat, "longitude")
  yy<-ncvar_get(dat, "latitude")
  varP <- ncvar_get(dat, "var__precipitation_rate__mm_hr_") 
  varSM <- ncvar_get(dat, "var__soil_moisture")
  nc_close(dat)
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whD<-which(as.Date(tms)==date-2)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  prec<-max(prec, max(varP[whx, why, whD]))
  whD<-which(as.Date(tms)==date-2)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  if(length(whD)==0) whD<-1
  soilm<-max(soilm, max(varSM[whx, why, whD]))
  return(c(sbr, prec, soilm))
}

get_prec_soilm<-function(date, lng, lat){
  sp <- SpatialPoints(cbind(lng, lat))
  doy <- as.numeric(strftime(date, format = "%j"))
  # 24 hours prior
  whD<-max(which(datesMET<=date-1))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
                                      collapse = ""),".nc")
  dat = nc_open(fln)
  tt  <- ncvar_get(dat, "time") 
  tms24<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  xx<-ncvar_get(dat, "longitude")
  yy<-ncvar_get(dat, "latitude")
  varP <- ncvar_get(dat, "var__precipitation_rate__mm_hr_") 
  varSM <- ncvar_get(dat, "var__soil_moisture")
  nc_close(dat)
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whD<-which(as.Date(tms24)==date-1)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  prec24<-varP[whx, why, whD]
  if(length(whD)==0) whD<-1
  soilm24<-varSM[whx, why, whD]
  tms24<-tms24[whD]
  # 48 hours prior
 whD<-max(which(datesMET<=date-2))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
                                      collapse = ""),".nc")
  dat = nc_open(fln)
  tt  <- ncvar_get(dat, "time") 
  tms48<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  xx<-ncvar_get(dat, "longitude")
  yy<-ncvar_get(dat, "latitude")
  varP <- ncvar_get(dat, "var__precipitation_rate__mm_hr_") 
  varSM <- ncvar_get(dat, "var__soil_moisture")
  nc_close(dat)
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whD<-which(as.Date(tms48)==date-2)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  prec48<-varP[whx, why, whD]
  if(length(whD)==0) whD<-1
  soilm48<-varSM[whx, why, whD]
  tms48<-tms48[whD]
  ans<-rbind(data.frame(time=tms48, value=prec48, type="Precipitation"),
             data.frame(time=tms48, value=soilm48, type="Soil moisture"),
             data.frame(time=tms24, value=prec24,type="Precipitation"),
             data.frame(time=tms24, value=soilm24, type="Soil moisture"))
  return(ans)
}

egg_development<-function(t){
  9.41*exp(-0.00357*(35.019-t)^2)*as.numeric(t>=10 & t<=34) 
}

egg_dev<-function(date, lng, lat){
  dev_day<-NA
  sp <- SpatialPoints(cbind(lng, lat))
  whD<-max(which(datesMET<=date))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
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
  dvlp_rt<-3*egg_development(ans$tempr)/24 # 3-hour development rate
  cm_dvlp_rt<-cumsum(dvlp_rt)
  wh<-which(cm_dvlp_rt>=100)
  if(length(wh)>0){
    wh<-min(wh)
    dev_day<-as.Date(ans$time[wh])
  } else{
    if(whD<length(datesMET)){
      fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+1]),"-")),
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
      dvlp_rt<-3*egg_development(ans$tempr)/24 # 3-hour development rate
      cm_dvlp_rt<-cumsum(dvlp_rt)
      wh<-which(cm_dvlp_rt>=100)
      if(length(wh)>0){
        wh<-min(wh)
        dev_day<-as.Date(ans$time[wh])
      } else {
        if(whD<length(datesMET)-1){
          fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+2]),"-")),
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
          dvlp_rt<-egg_development(ans$tempr)*3/24 # 3-hour temperature data
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
  if (as.numeric(dev_day-date)<10) dev_day<-NA
  if (as.numeric(dev_day-date)>65) dev_day<-NA
  return(dev_day)
}

temp_profile_egg_dev<-function(date, lng, lat){
  dev_day<-NA
  sp <- SpatialPoints(cbind(lng, lat))
  whD<-max(which(datesMET<=date))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
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
  dvlp_rt<-3*egg_development(ans$tempr)/24 # 3-hour development rate
  cm_dvlp_rt<-cumsum(dvlp_rt)
  wh<-which(cm_dvlp_rt>=100)
  if(length(wh)>0){
    wh<-min(wh)
    dev_day<-as.Date(ans$time[wh])
  } else{
    if(whD<length(datesMET)){
      fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+1]),"-")),
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
      dvlp_rt<-3*egg_development(ans$tempr)/24 # 3-hour development rate
      cm_dvlp_rt<-cumsum(dvlp_rt)
      wh<-which(cm_dvlp_rt>=100)
      if(length(wh)>0){
        wh<-min(wh)
        dev_day<-as.Date(ans$time[wh])
      } else {
        if(whD<length(datesMET)-1){
          fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+2]),"-")),
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
          dvlp_rt<-egg_development(ans$tempr)*3/24 # 3-hour temperature data
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
  colnames(ans)<-c("time", "value")
  return(ans)
}

hopper_development<-function(t){
  (0.222*t -3.166)*as.numeric(t>=23 & t<=32) 
}


smooth_Whittaker<-function(y, lambda=10^1.3){
  m = length(y)
  w<-rep(1, m)
  D = diff(diag.spam(m), diff = 2)
  W = diag.spam(w)
  z = solve(W + lambda * t(D) %*% D, w * y)
  return(z)
}

fit_line<-function(date1, date2, val1, val2){
  x<-c(0, as.numeric(date2-date1))
  y<-c(val1, val2)
  slope <- diff(y)/diff(x)
  intercept <- y[1]-slope*x[1]
  return(data.frame( slope,  intercept))
}


find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


get_NDVI_values<-function(lng, lat, day1, day2){
  sp <- SpatialPoints(cbind(lng, lat))
  pref<-"data/NDVI/AH_NDVI_16d_1km_"
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
  ans<-ans[order(ans$date),]
  ans<-ans[which(ans$date>=day1 & ans$date<=day2),]
  return(ans)
}

hopper_dev_dur<-function(date, lng, lat){
  dev_day<-NA
  sp <- SpatialPoints(cbind(lng, lat))
  whD<-max(which(datesMET<=date))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
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
  dvlp_rt<-hopper_development(ans$tempr)*3/24 # 3-hour intervals
  cm_dvlp_rt<-cumsum(dvlp_rt)
  wh<-which(cm_dvlp_rt>=100)
  if(length(wh)>0){
    wh<-min(wh)
    dev_day<-as.Date(ans$time[wh])
  } else{
    if(whD<length(datesMET)){
      fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+1]),"-")),
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
      ans<-rbind(ans, data.frame(time=tms[whDD], tempr=varT[whx, why,whDD]))
      dvlp_rt<-hopper_development(ans$tempr)*3/24 # 3-hour intervals
      cm_dvlp_rt<-cumsum(dvlp_rt)
      wh<-which(cm_dvlp_rt>=100)
      if(length(wh)>0){
        wh<-min(wh)
        dev_day<-as.Date(ans$time[wh])
      } else {
        if(whD<length(datesMET)-1){
          fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+2]),"-")),
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
  if (as.numeric(dev_day-date)<24) dev_day<-NA
  if (as.numeric(dev_day-date)>95) dev_day<-NA
  return(dev_day)
}

temp_profile_hopper_dev_dur<-function(date, lng, lat){
  dev_day<-NA
  sp <- SpatialPoints(cbind(lng, lat))
  whD<-max(which(datesMET<=date))
  fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD]),"-")),
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
  dvlp_rt<-hopper_development(ans$tempr)*3/24 # 3-hour intervals
  cm_dvlp_rt<-cumsum(dvlp_rt)
  wh<-which(cm_dvlp_rt>=100)
  if(length(wh)>0){
    wh<-min(wh)
    dev_day<-as.Date(ans$time[wh])
  } else{
    if(whD<length(datesMET)){
      fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+1]),"-")),
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
      ans<-rbind(ans, data.frame(time=tms[whDD], tempr=varT[whx, why,whDD]))
      dvlp_rt<-hopper_development(ans$tempr)*3/24 # 3-hour intervals
      cm_dvlp_rt<-cumsum(dvlp_rt)
      wh<-which(cm_dvlp_rt>=100)
      if(length(wh)>0){
        wh<-min(wh)
        dev_day<-as.Date(ans$time[wh])
      } else {
        if(whD<length(datesMET)-1){
          fln<-paste0("data/Met_Data_", str_c(unlist(strsplit(as.character(datesMET[whD+2]),"-")),
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
  colnames(ans)<-c("time", "value")
  return(ans)
}

hopper_dev<-function(date, lng, lat, thresh, peak=TRUE, lambda=1, m=2){
  ans<-NA
  dev.day<-hopper_dev_dur(date, lng, lat)
  if(!is.na(dev.day)){
    dev.per<-as.numeric(dev.day-date)
    if(dev.per>=24 & dev.per<=95){
      day1<-date
      day2<-dev.day
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
  }
  return(ans)
}

get_wind<-function(date, lng, lat, rn.sd){
  tmp<-sqrt((src$lng-lng)^2 +(src$lat-lat)^2)
  wh<-which(tmp==min(tmp))
  d.lng<-src$lng[wh]-lng
  d.lat<-src$lat[wh]-lat
  if(wh<10) s<-paste0("0000",wh)
  if(wh>=10 & wh<100) s<-paste0("000",wh)
  if(wh>=100 & wh<1000) s<-paste0("00",wh)
  if(wh>=1000 & wh<10000) s<-paste0("0",wh)
  dt<-format(date,"%Y%m%d")
  flnm <- paste0("data/wind_trajectories/Data_Traj_Swarm", s, "_C1_",dt,".nc")
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

get_wind_ind<-function(date, lng, lat, ind){
  tmp<-sqrt((src$lng-lng)^2 +(src$lat-lat)^2)
  wh<-which(tmp==min(tmp))
  d.lng<-src$lng[wh]-lng
  d.lat<-src$lat[wh]-lat
  if(wh<10) s<-paste0("0000",wh)
  if(wh>=10 & wh<100) s<-paste0("000",wh)
  if(wh>=100 & wh<1000) s<-paste0("00",wh)
  if(wh>=1000 & wh<10000) s<-paste0("0",wh)
  dt<-format(date,"%Y%m%d")
  flnm <- paste0("data/wind_trajectories/Data_Traj_Swarm", s, "_C1_",dt,".nc")
  dat = nc_open(flnm)
  tt <- ncvar_get(dat, "time")            # Extract time scale
  tms<-utcal.nc("hours since 1970-01-01 00:00:00", tt, "c")
  yy<-ncvar_get(dat, "latitude")
  xx<-ncvar_get(dat, "longitude")
  nc_close(dat)
  j<-ind
  wn.trj<-data.frame(date=tms, lng=xx[j,]-d.lng, lat=yy[j,]-d.lat)
  wn.trj<-wn.trj[order(wn.trj$date),]
  return(wn.trj)
}

get_land_cover<-function(lng, lat){
  lncv <- read.csv("data/Land_cover_1km.csv", header = TRUE)
  ds<-sqrt((lncv$lng-lng)^2+(lncv$lat-lat)^2)
  wh<-which(ds==min(ds))
  ans<-lncv$Type[wh]
  rm(lncv)
  return(ans)
}

get_NDVI_trend<-function(yy, npsi=6){ # xx<-ndvi data
  xx<-yy
  day1<-min(xx$date); day2<-max(xx$date)
  for(i in 2:nrow(xx)){
    cfs<-fit_line(xx$date[i-1], xx$date[i], xx$NDVIsmooth[i-1],
                  xx$NDVIsmooth[i])
    ds<-seq(as.Date(xx$date[i-1]),as.Date(xx$date[i]), by=1)
    for(j in 2:(length(ds)-1)){
      val<-cfs$intercept+cfs$slope*as.numeric(ds[j]-ds[1])
     xx<-rbind(xx, data.frame(date=as.Date(ds[j]), NDVI=1,
                              NDVIsmooth=val))
    }
  }
  xx<-xx[order(xx$date),]
  xx<-xx[which(xx$date>=day1 & xx$date<=day2),]
  ans<-"Constant"
  xx<-xx[rev(order(xx$date)),]
  x<-data.frame(x=1:nrow(xx), y=xx$NDVIsmooth)
  fit_lm = lm(y ~ 1 + x, data = x)  # intercept-only model
  fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi =  npsi)  # Two change points along x
  ch.pnt<-round(as.numeric(fit_segmented$psi[,2]))
  if(length(ch.pnt)==0){
    fit_segmented = segmented(fit_lm, seg.Z = ~x, npsi =  2)  # Two change points along x
    ch.pnt<-round(as.numeric(fit_segmented$psi[,2]))
  }
  tl<-ch.pnt[1]
  x0<-x$y[1]
  x1<-x$y[tl]
  if((x0<x1) & (abs(x0-x1)>0.1)) ans<-"Decreasing"
  if((x1<x0) & (abs(x0-x1)>0.1)) ans<-"Increasing"
  if(length(ch.pnt)>1){
    x2<-x$y[ch.pnt[2]]
    if(ans=="Decreasing" & x1<x2) {
      x1<-x2; 
      tl<-ch.pnt[2]
      if(length(ch.pnt)>2){
        x3<-x$y[ch.pnt[3]]
        if(x2<x3) {x1<-x3; tl<-ch.pnt[3]}
      }}
    if(ans=="Increasing" & x2<x1) {
      x1<-x2
      tl<-ch.pnt[2]
      if(length(ch.pnt)>2){
        x3<-x$y[ch.pnt[3]]
        if(x3<x2) {x1<-x3; tl<-ch.pnt[3]}
      }}
  }
  
  return(data.frame(trend=ans, days=tl, x_start=x1, x_end=x0))
}

get_NDVI_dnst<-function(x){
  ans<- "No vegetation" #(NDVI<0), 
  if(x>0 & x<=0.15) ans<- "Lowest density" # 
  if(x>0.15 & x<=0.3) ans<- "Lower density" #
  if(x>0.3 & x<=0.45) ans<- "Dense vegetation"
  if(x>0.45 & x<=0.6) ans<- "Higher density"
  if(x>0.6) ans<-"Highest density"
  return(ans)
}


get_lct<-function(lct){
  ans<-"Other"
  if(lct %in% c(20, 30, 90)) ans<-"Shrubs" # Shrubs; Herbaceous vegetation; Herbaceous wetland
  if(lct==40) ans<-"Cropland" # Cultivated and managed vegetation/agriculture
  if(lct==60) ans<-"Sparse vegetation" # Bare / sparse vegetation
  if(lct==50) ans<-"Other" # Urban/built up
  if(lct %in% c(111,112, 113, 114, 115, 116, 121, 122, 123,
                124, 125, 126)) ans<-"Forest" #  Forest
  if(lct==70) ans<-"Other" # Snow and ice
  if(lct %in% c(0, 80, 200)) ans<-"Other" # Unknown, Permanent water bodies or Oceans, seas
  if(lct==100) ans<-"Sparse vegetation" # Moss and lichen
  return(ans)
}


# Land cover types
vgts.lct<-c("Sparse vegetation", "Shrubs", "Forest", "Cropland")
vgts.dnst<-c("Lowest density", "Lower density",
             "Dense vegetation", "Higher density", "Highest density")
vgts.trnd<-c("Decreasing", "Constant", "Increasing" )
st.lngth<-rbind(data.frame(length="Short stay",  days.min=1, days.max=2),
                data.frame(length="Medium stay",  days.min=2, days.max=4),
                data.frame(length="Long stay",  days.min=4, days.max=7))

lng.stay<-data.frame(lc=c(), tr=c(), vl=c(), stay=c())
for(v1 in vgts.lct){
  for(v2 in vgts.trnd){
    for(v3 in vgts.dnst){
      lng.stay<-rbind(lng.stay,
                      data.frame(lc=v1, tr=v2, vl=v3, stay="Short stay"))
    }
  }
}

lng.stay$stay<-NA
# Long stay
lng.st<-rbind(data.frame(lc="Cropland", tr=c("Decreasing", "Constant", "Increasing" ), vl="Highest density"),
              data.frame(lc="Cropland", tr="Constant",vl="Higher density"),
              data.frame(lc="Cropland", vl="Higher density", tr="Increasing"),
              data.frame(lc="Forest", vl="Higher density", tr="Increasing"),
              data.frame(lc="Forest", vl="Highest density", tr="Increasing"))
for(ii in 1:nrow(lng.st)){
  wh<-which(lng.stay$lc==lng.st$lc[ii] & lng.stay$tr==lng.st$tr[ii] & lng.stay$vl==lng.st$vl[ii])
  if(length(wh)>0) lng.stay$stay[wh]<-"Long stay"
}
#  Medium stay
mdm.st<-rbind(data.frame(lc="Cropland", vl=c("Lower density", "Dense vegetation"), tr="Increasing"),
              data.frame(lc="Cropland", vl=c("Dense vegetation", "Higher density"), tr="Decreasing"),
              data.frame(lc="Cropland", vl="Dense vegetation", tr="Constant"),
              data.frame(lc="Shrubs", vl="Lower density", tr="Increasing"),
              data.frame(lc="Shrubs", vl="Dense vegetation", tr=c("Increasing","Constant")),
              data.frame(lc="Forest", vl="Lower density", tr="Increasing"),
              data.frame(lc="Forest", vl="Dense vegetation", tr=c("Increasing","Constant")),
              data.frame(lc="Forest", vl="Higher density", tr=c("Constant","Decreasing")),
              data.frame(lc="Sparse vegetation", vl="Dense vegetation", tr="Increasing"),
              data.frame(lc="Forest", vl="Highest density", tr=c("Constant","Decreasing")),
              data.frame(lc="Shrubs", vl="Higher density", tr=c("Increasing","Constant")))
for(ii in 1:nrow(mdm.st)){
  wh<-which(lng.stay$lc==mdm.st$lc[ii] & lng.stay$tr==mdm.st$tr[ii] & lng.stay$vl==mdm.st$vl[ii])
  if(length(wh)>0) lng.stay$stay[wh]<-"Medium stay"
}
#  Short stay
shr.st<-rbind(data.frame(lc="Cropland", vl=c( "Lower density"), tr= "Constant"),
              data.frame(lc="Cropland", vl=c( "Lowest density", "Lower density"), tr= "Constant"),
              data.frame(lc="Cropland", vl="Lowest density", tr="Increasing"),
              data.frame(lc="Cropland", vl="Lower density", tr="Constant"),
              data.frame(lc="Cropland", vl="Lower density", tr="Decreasing"),
              data.frame(lc="Shrubs", vl="Lowest density", tr=c("Increasing","Constant")),
              data.frame(lc="Shrubs", vl="Lower density", tr=c("Decreasing","Constant")),
              data.frame(lc="Shrubs", vl="Dense vegetation", tr="Decreasing"),
              data.frame(lc="Sparse vegetation", vl="Lowest density",tr=c("Increasing","Constant","Decreasing")),
              data.frame(lc="Sparse vegetation", vl="Lower density",tr=c("Increasing","Constant","Decreasing")),
              data.frame(lc="Forest", vl="Lowest density", tr=c("Increasing","Constant")),
              data.frame(lc="Forest", vl="Lower density", tr=c("Decreasing","Constant")),
              data.frame(lc="Forest", vl="Dense vegetation", tr="Decreasing"),
              data.frame(lc=c("Cropland", "Shrubs","Forest"), vl="Lowest density", tr="Decreasing"),
              data.frame(lc="Sparse vegetation", vl="Dense vegetation", tr=c("Decreasing","Constant")),
              data.frame(lc="Shrubs", vl="Higher density", tr="Decreasing"))
for(ii in 1:nrow(shr.st)){
  wh<-which(lng.stay$lc==shr.st$lc[ii] & lng.stay$tr==shr.st$tr[ii] & lng.stay$vl==shr.st$vl[ii])
  if(length(wh)>0) lng.stay$stay[wh]<-"Short stay"
}

lng.stay$lc<-factor(lng.stay$lc, levels = vgts.lct)
lng.stay$tr<-factor(lng.stay$tr, levels = rev(vgts.trnd))
lng.stay$vl<-factor(lng.stay$vl, levels = rev(vgts.dnst))
lng.stay$stay<-factor(lng.stay$stay, levels = st.lngth$length)

lng.stay$VL<-"TBA"
lng.stay$VL[lng.stay$vl=="Lowest density"]<-"Lowest density: (0,0.15]"
lng.stay$VL[lng.stay$vl=="Lower density"]<-"Lower density: (0.15,0.3]"
lng.stay$VL[lng.stay$vl=="Dense vegetation"]<-"Dense vegetation: (0.3,0.45]"
lng.stay$VL[lng.stay$vl=="Higher density"]<-"Higher density: (0.45,0.6]"
lng.stay$VL[lng.stay$vl=="Highest density"]<-"Highest density: (0.6,1]"


swrm_stay<-function(lct, ndvi, npsi, rnm.sd=c()){
  ans<-1
  vi.lct<-get_lct(lct)
  if(!vi.lct=="Other"){
    vi.tr<-get_NDVI_trend(ndvi, npsi)
    vi.dst<-get_NDVI_dnst(vi.tr$x_end)
    wh<-which(lng.stay$lc==vi.lct & 
                lng.stay$tr==vi.tr$trend & 
                lng.stay$vl== vi.dst)
    if(length(wh)==1){
      st.nm<-lng.stay$stay[wh]
      if(!is.na(st.nm)){
        wh<-which(st.lngth$length==st.nm)
        if(length(rnm.sd)>0) set.seed(rnm.sd)
        ans<-sample(seq(st.lngth$days.min[wh], 
                        st.lngth$days.max[wh], by=1), 1)
      }
    }
  }
  return(ans)
}


