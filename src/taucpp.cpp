#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List taucpp(Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector xstatus, 
            Rcpp::NumericVector ystatus, Rcpp::NumericVector KMxsurv, 
            Rcpp::NumericVector KMysurv, Rcpp::NumericVector KMxtime, 
            Rcpp::NumericVector KMytime){
  int n = x.size();
  Rcpp::NumericMatrix a(n);
  Rcpp::NumericMatrix b(n);
  for(int i=0;i<n;++i){
    for(int j=0;j<i;++j){
      if((xstatus(i) == 1) & (xstatus(j) == 1)){
        a(i,j) = (x(i) > x(j)) - (x(i) < x(j));
      }
      else if((xstatus(i) == 0) & (xstatus(j) == 1)){
        Rcpp::NumericVector helperj = KMxsurv[KMxtime == x(j)];
        Rcpp::NumericVector helperi = KMxsurv[KMxtime == x(i)];
        a(i,j) = (x(i) >= x(j)) + (x(i) < x(j)) *
	        (2 * sum(helperj) / sum(helperi) - 1);
      }
      else if((xstatus(i) == 1) & (xstatus(j) == 0)){
        Rcpp::NumericVector helperj = KMxsurv[KMxtime == x(j)];
        Rcpp::NumericVector helperi = KMxsurv[KMxtime == x(i)];
        a(i,j) = (x(i) > x(j)) * (1 - 2 * sum(helperi) / sum(helperj)) - (x(i) <= x(j));
      }
      else{
        Rcpp::NumericVector helperj = KMxsurv[KMxtime == x(j)];
        Rcpp::NumericVector helperi = KMxsurv[KMxtime == x(i)];
        a(i,j) = (x(i) > x(j)) * (1 - sum(helperi) / sum(helperj)) + (x(i) < x(j)) * 
	        (sum(helperi) / sum(helperj) - 1);
      }

      if((ystatus(i) == 1) & (ystatus(j) == 1)){
        b(i,j) = (y(i) > y(j)) - (y(i) < y(j));
      }
      else if((ystatus(i) == 0) & (ystatus(j) == 1)){
        Rcpp::NumericVector helperj = KMysurv[KMytime == y(j)];
        Rcpp::NumericVector helperi = KMysurv[KMytime == y(i)];
        b(i,j) = (y(i) >= y(j)) + (y(i) < y(j)) * (2 * sum(helperj) / sum(helperi) - 1);
      }
      else if((ystatus(i) == 1) & (ystatus(j) == 0)){
        Rcpp::NumericVector helperj = KMysurv[KMytime == y(j)];
        Rcpp::NumericVector helperi = KMysurv[KMytime == y(i)];
        b(i,j) = (y(i) > y(j)) * (1 - 2 * sum(helperi) / sum(helperj)) - (y(i) <= y(j));
      }
      else{
        Rcpp::NumericVector helperj = KMysurv[KMytime == y(j)];
        Rcpp::NumericVector helperi = KMysurv[KMytime == y(i)];
        b(i,j) = (y(i) > y(j)) * (1 - sum(helperi) / sum(helperj)) + (y(i) < y(j)) * 
          (sum(helperi) / sum(helperj) - 1);
      }
    }
  }
  Rcpp::NumericVector c = a * b;
  return Rcpp::List::create(Rcpp::Named("tau") = sum(c) / pow(sum(pow(a, 2)) * 
                      sum(pow(b, 2)), .5), Rcpp::Named("a") = a, Rcpp::Named("b") = b);
}

