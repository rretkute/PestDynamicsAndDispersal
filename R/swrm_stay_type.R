#' Get type of stay based on land cover type and NDVI profile.
#'
#' @param lct Land cover type.
#' @param ndvi  Vector of NDVI values.
#' @param npsi Number of change points in NDVI profile.
#' @return Type of stay: short, medium or long.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

swrm_stay_type<-function(lct, ndvi, npsi=8){
  ans<-"tba"
  vi.lct<-get_lct(lct)
  if(!vi.lct=="Other"){
    vi.tr<-get_NDVI_trend(ndvi, npsi)
    vi.dst<-get_NDVI_dnst(vi.tr$x_end)
    wh<-which(lng.stay$lc==vi.lct & 
                lng.stay$tr==vi.tr$trend & 
                lng.stay$vl== vi.dst)
    ans<-as.character(lng.stay$stay[wh])
  }
  return(ans)
}