#' Length of time a swarm stays at a location.
#' 
#'
#' @param lct Land cover type.
#' @param ndvi Profile of NDVI.
#' @param npsi Number of change points.
#' @return Number of days.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

swrm_stay<-function(lct, ndvi, npsi){
  ans<-0
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
        ans<-sample(seq(st.lngth$days.min[wh], 
                        st.lngth$days.max[wh], by=1), 1)
      }
    }
  }
  return(ans)
}