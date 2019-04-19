#include <Rcpp.h>
using namespace Rcpp;

//' Matrix with expected values for CHR independence test
//'
//' @title Matrix with expected values for CHR independence test
//' @param x,y Vectors of failure times
//' @param xstatus,ystatus Status indicators for failure times
//' @param xuni,yuni Grid points to get expected value for
//' @return Matrix with expected values for CHR independence test
//' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com
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
