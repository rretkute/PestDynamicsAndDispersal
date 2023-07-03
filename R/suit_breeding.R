#' Calcualte suitability for breeding/laying eggs
#'
#' @param date Date
#' @param lng Longitude
#' @param lat Latitude
#' @param sbr Raster with breeding suitability scores
#' @param METdata  MET Office data on soil moisture and precipitation
#' @return  Values of breeding suitability, total precipitation and max soil moisture
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
suit_breeding<-function(date, lng, lat, sbr, METdata){
  xx<-METdata[[1]]
  yy<-METdata[[2]]
  tms<-METdata[[3]]
  varP<-METdata[[4]]
  varSM<-METdata[[5]]
  METdata<-list(xx, yy, tms, varP, varSM)
  sp <- SpatialPoints(cbind(lng, lat))
  doy <- as.numeric(strftime(date, format = "%j"))
  # Breeding suitability map
  sbr<-raster::extract(sbr, sp, method='bilinear')/100
  if(is.na(sbr)) sbr<-0
  whx<-which((xx-lng)^2==min((xx-lng)^2))[1]
  why<-which((yy-lat)^2==min((yy-lat)^2))[1]
  whD<-which(as.Date(tms)==date)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  prec<-3*sum(varP[whx, why, whD])
  whD<-which(as.Date(tms)==date)  # @ 2, 5, 8, 11, 14, 17, 20 & 23
  if(length(whD)==0) whD<-1
  soilm<-max(varSM[whx, why, whD])
  return(c(sbr, prec, soilm))
}