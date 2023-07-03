#' Get NDVI trend.
#'
#' @param xx NDVI dates and values.
#' @param ndvi  Vector of NDVI values.
#' @param npsi Number of change points in NDVI profile.
#' @return Trend, change date and values at the change pint and at the day.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

get_NDVI_trend<-function(xx, npsi=10){ # ndvi data
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
  if(x0<x1) ans<-"Decreasing"
  if(x1<x0) ans<-"Increasing"
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