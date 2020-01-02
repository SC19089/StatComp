#include <Rcpp.h>
using namespace Rcpp;

//' @title etimate possion process.
//' @description we can get in in R.
//' @param lambda the parameter of the exponential distribution.
//' @return a result \code{n}
//' @examples
//' \dontrun{
//' xtC(100,1)
//' }
//' @export
// [[Rcpp::export]]
NumericVector xtC(int N,double lambda){
  NumericVector x=rexp(N,lambda);
  NumericVector Tn(N);
  NumericVector xt(1);
  Tn[0]=x[0];
  for(int i=1;i<N;i++){
    Tn[i]=Tn[i-1]+x[i];
  } 
  for(int i=0;i<N;i++){
    int t=10;
    if(Tn[i]<=t&&Tn[i+1]>t){
      xt=i;}}
  return xt;
}