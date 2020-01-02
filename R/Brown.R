
#' @title Brownian Motion.
#' @description plot of st Brownian Motion，before use the function,you should use package e1071.
#' @param t the time.
#' @param n the number of paces of Brownian Motion.
#' @return plot figure\code{n}
#' @examples
#' \dontrun{
#' Brown(1,6)
#' }
#' @export
Brown=function(t,n){
  N=1000
  m=n-1
  time=seq(0,t,length=N)
  path=e1071::rwiener(end=1,frequency=N)
  plot(time,path,xlab="time",type="l",main=" Standard BM")
  for(i in 1:m){
    path=e1071::rwiener(end=1,frequency=N)
    lines(time,path)
  }
}


