#' Sample hopper development period from the truncated lognormal distribution
#'
#' @param meanlog Mean length of period (default is 36 days).
#' @param sdlog Standart deviation of period (default is q days).
#' @return  Length of egg development period in days
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
sample_dev_period<-function(meanlog = log(36), sdlog =1){
  cnd<-TRUE
  while(cnd){
    nn<-round(rlnorm(1, meanlog = meanlog, sdlog =sdlog))+24
    if(nn>=24 & nn<=96) cnd<-FALSE
  }
  return(nn)
}