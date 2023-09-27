#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube gghelper(NumericMatrix Lambda, NumericVector time, 
              NumericVector birthYear, NumericMatrix R, NumericMatrix xi, 
              NumericMatrix Y, double theta, double beta, double n0, 
              double n1){
  cube g = cube(n0, n0, n1);
  for(int u = 0; u < n0; ++u){
    std::cout << u << std::endl;
    for(int t = 0; t < n0;++t){
      for(int i = 0; i < n1; ++i){
        g(u, t, i) = theta / R(t,i) * xi(t,i) * 
          (exp(theta * Lambda(t,2*i) * exp(beta * birthYear(2*i))) *
          exp(beta * birthYear(2*i)) * Y(u,2*i) +
          exp(theta * Lambda(t,2*i+1) * exp(beta * birthYear(2*i+1))) *
          exp(beta * birthYear(2*i+1)) * Y(u,2*i+1));
      }
    }
  }
  return g;
}
