#' Fit a line of NDVI values between two dates.
#'
#' @param date1 Date 1.
#' @param date2 Date 2.
#' @param val1  NDVI value at date1.
#' @param val2 NDVI value at date2.
#' @return Smoothed values.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

fit_line<-function(date1, date2, val1, val2){
  x<-c(0, as.numeric(date2-date1))
  y<-c(val1, val2)
  slope <- diff(y)/diff(x)
  intercept <- y[1]-slope*x[1]
  return(data.frame( slope,  intercept))
}
