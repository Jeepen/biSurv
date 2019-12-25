#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix martin(NumericVector time, NumericVector status, NumericVector times,
		     NumericVector cumhaz, int n, int nu){
  NumericMatrix M(nu,n);
  for(int t=0;t<nu;++t){
    for(int i=0;i<n;++i){
      NumericVector help = cumhaz[times == time(i)];
      M(t,i) = status(i) * (time(i) <= times(t)) - cumhaz(t) * (time(i) > times(t)) -
	sum(help) * (time(i) <= times(t));
    }
  }
  return M;
} 
