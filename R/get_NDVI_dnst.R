#' Get vegetation density type.
#'
#' @param xx NDVI value.
#' @return Density type (Lowest, Lower, Dense, Higher or Highest).
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

get_NDVI_dnst<-function(x){
  ans<- "No vegetation" #(NDVI<0), 
  if(x>0 & x<=0.15) ans<- "Lowest density" # 
  if(x>0.15 & x<=0.3) ans<- "Lower density" #
  if(x>0.3 & x<=0.45) ans<- "Dense vegetation"
  if(x>0.45 & x<=0.6) ans<- "Higher density"
  if(x>0.6) ans<-"Highest density"
  return(ans)
}