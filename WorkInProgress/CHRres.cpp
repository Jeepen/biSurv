#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericMatrix CHRresCpp(NumericVector x, NumericVector y, NumericVector xstatus,
			NumericVector ystatus, NumericVector xuni, NumericVector yuni){
  double n1 = xuni.size();
  double n2 = yuni.size();
  arma::mat R(n1, n2); 
  arma::mat L11(n1,n2); 
  arma::mat L10(n1,n2);
  arma::mat L01(n1,n2);
  for(int i=0; i<n1; ++i){
    for(int j=0; j<n2; ++j){
      R(i,j) = sum((x >= xuni(i)) % (y >= yuni(j)));
      L11(i,j) = sum(xstatus % ystatus % ((x == xuni(i)) % (y == yuni(j)))) /
        R(i,j);
      S(i,j) = xstatus;
    }
  }
}
