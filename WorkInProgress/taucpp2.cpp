#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List taucpp2(const arma::vec &x, const arma::vec &y, const arma::vec &xstatus, 
            const arma::vec &ystatus, const arma::vec &KMxsurv, 
		   const arma::vec &KMysurv, const arma::vec &KMxtime, 
            const arma::vec &KMytime){
  int n = x.size();
  // arma::mat a(n,n);
  // arma::mat b(n,n);
  Rcpp::NumericMatrix a(n);
  Rcpp::NumericMatrix b(n);
  for(int i=0;i<n-1;++i){
    for(int j=0;j<i-1;++j){
      if((xstatus(i) == 1) & (xstatus(j) == 1)){
        a(i,j) = (x(i) > x(j)) - (x(i) < x(j));
      }
      else if((xstatus(i) == 0) & (xstatus(j) == 1)){
        a(i,j) = (x(i) >= x(j)) + (x(i) < x(j)) *
	  (2 * as_scalar(KMxsurv(find(KMxtime == x(j)))) / 
	   as_scalar(KMxsurv(find(KMxtime == x(i)))) - 1);
      }
      else if((xstatus(i) == 1) & (xstatus(j) == 0)){
        a(i,j) = (x(i) > x(j)) * (1 - 2 * as_scalar(KMxsurv(find(KMxtime == x(i)))) / 
				  as_scalar(KMxsurv(find(KMxtime == x(j))))) - (x(i) <= x(j));
      }
      else{
        a(i,j) = (x(i) > x(j)) * (1 - as_scalar(KMxsurv(find(KMxtime == x(i)))) / 
				  as_scalar(KMxsurv(find(KMxtime == x(j))))) + (x(i) < x(j)) * 
	        (as_scalar(KMxsurv(find(KMxtime == x(i)))) / 
		 as_scalar(KMxsurv(find(KMxtime == x(j)))) - 1);
      }

      if((ystatus(i) == 1) & (ystatus(j) == 1)){
        b(i,j) = (y(i) > y(j)) - (y(i) < y(j));
      }
      else if((ystatus(i) == 0) & (ystatus(j) == 1)){
        b(i,j) = (y(i) >= y(j)) + (y(i) < y(j)) * 
	  (2 * as_scalar(KMysurv(find(KMytime == y(j)))) / 
	   as_scalar(KMysurv(find(KMytime == y(i)))) - 1);
      }
      else if((ystatus(i) == 1) & (ystatus(j) == 0)){
        b(i,j) = (y(i) > y(j)) * (1 - 2 * as_scalar(KMysurv(find(KMytime == y(i)))) / 
				  as_scalar(KMysurv(find(KMytime == y(j))))) - (y(i) <= y(j));
      }
      else{
        b(i,j) = (y(i) > y(j)) * (1 - as_scalar(KMysurv(find(KMytime == y(i)))) / 
				  as_scalar(KMysurv(find(KMytime == y(j))))) + (y(i) < y(j)) * 
          (as_scalar(KMysurv(find(KMytime == y(i)))) / as_scalar(KMysurv(find(KMytime == y(j)))) - 1);
      }
    }
  }
  Rcpp::NumericVector c = a * b;
  // arma::mat c = a % b;
  return Rcpp::List::create(Rcpp::Named("tau") = sum(c) / pow(sum(pow(a, 2)) * 
                      sum(pow(b, 2)), .5), Rcpp::Named("a") = a, Rcpp::Named("b") = b);
}

