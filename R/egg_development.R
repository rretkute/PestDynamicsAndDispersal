#' Egg development rate  as a function of temperature
#' .
#'
#' @param t Average daily temperature.
#' @return Development rate (per day).
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

egg_development<-function(t){
  9.41*exp(-0.00357*(35.019-t)^2)*as.numeric(t>10)
}