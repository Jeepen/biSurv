#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]] 
NumericMatrix eyyfunc(NumericVector x, NumericVector y, NumericVector xuni, 
                      NumericVector yuni){
  int nx = xuni.size();
  int ny = yuni.size();
  NumericMatrix eyy(nx,ny);
  for(int i=0; i<nx; ++i){
    for(int j=0; j<ny; ++j){
      eyy(i,j) = mean((x >= xuni(i)) & (y >= yuni(j)));
    }
  }
  return eyy;
}
