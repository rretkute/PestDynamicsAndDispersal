#' Whittaker smoothing for NDVI values
#'
#' @param y Vector of NDVI values.
#' @param lambda Smoothing window.
#' @return Smoothed values.
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

smooth_Whittaker<-function(y, lambda=10^1.3){
  m = length(y)
  w<-rep(1, m)
  D = diff(diag.spam(m), diff = 2)
  W = diag.spam(w)
  z = solve(W + lambda * t(D) %*% D, w * y)
  return(z)
}