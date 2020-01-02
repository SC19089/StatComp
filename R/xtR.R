
#' @title Creat a possion process 
#' @description Estimate the number of the events when t=10.
#' @param N the number of Tn.
#' @param lambda the parameter of the exponential distribution.
#' @return a estimate result \code{n}
#' @examples
#' \dontrun{
#' re2R=numeric(1e4)
#' for(i in 1:1e4){
#' re2R[i] = xtR(100,1)}
#' round(mean(re2R))
#' }
#' @export
xtR=function(N,lambda){
  x=rexp(N,lambda)
  Tn=vector(length=N)
  Tn[1]=x[1]
  for(i in 2:N){
    Tn[i]=Tn[i-1]+x[i]
  } #the total time
  t=10
  i=1
  for(i in 1:N){
    if(Tn[i]<=t&&Tn[i+1]>t){
      xt=i
    }
  }
  xt
}


