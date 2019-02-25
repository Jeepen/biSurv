#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List hazardscpp(arma::vec x, arma::vec y, arma::vec xstatus, 
                  arma::vec ystatus, arma::vec xuni, 
                  arma::vec yuni){
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
      L10(i,j) = sum(xstatus % ((x == xuni(i)) % (y >= yuni(j)))) / R(i,j);
      L01(i,j) = sum(ystatus % ((x >= xuni(i)) % (y == yuni(j)))) / R(i,j);
    }
  }
  return Rcpp::List::create(Rcpp::Named("lambda11") = L11, 
                            Rcpp::Named("lambda10") = L10, 
                            Rcpp::Named("lambda01") = L01);
}

