
#' @title Use stratified sampling to estimate the integration.
#' @description Estimate the integration of x^nexp(-x) from 0 to 1.
#' @param n the power of x.
#' @return estimate value \code{n}
#' @examples
#' \dontrun{
#' re1R = sseR(5)
#' }
#' @export
sseR=function(n){
  f=function(x){x^n*exp(-x)*(x>0)*(x<1)}
  s=numeric(4)
  for(i in 1:4){
    s[i]=mean(f(runif(4,(i-1)*0.25,i*0.25)))
  }
  mean(s)
}




