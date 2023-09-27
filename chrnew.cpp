#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List chrResCpp(NumericVector x, NumericVector y, NumericVector xstatus, NumericVector ystatus, 
	       NumericVector xuni, NumericVector yuni, NumericMatrix W){
  // Initialize
  int n1 = xuni.size();
  int n2 = yuni.size();
  int nx = x.size();
  int ny = y.size();
  NumericVector lambdax(n1);
  NumericVector lambday(n2);
  NumericMatrix dNx(n1, nx);
  NumericMatrix Yx(n1, nx);
  NumericMatrix dMx(n1,nx);
  NumericMatrix dNy(n2, ny);
  NumericMatrix Yy(n2, ny);
  NumericMatrix dMy(n2,ny);
  NumericMatrix R(n1,n2);
  // R
  /*
  for(int i=0;i<n1;++i){
    for(int j=0;j<n2;++j){
      R(i,j) = sum(x >= xuni(i) & y >= yuni(j));
    }
  }
  */
  // lambdax
  for(int i=0; i<n1;++i){
    lambdax(i) = sum(x == xuni(i) & xstatus == 1) / sum(x>=xuni(i));
  }
  // lambday
  for(int i=0; i<n2;++i){
    lambday(i) = sum(y == yuni(i) & ystatus == 1) / sum(y>=yuni(i));
  }
  // dNx & Yx
  for(int i=0;i<n1;++i){
    for(int j=0;j<nx;++j){
      dNx(i,j) = (x(j) == xuni(i) & xstatus(j) == 1);
      Yx(i,j) = (x(j) >= xuni(i));
    }
  }
  // dNy & Yy
  for(int i=0;i<n2;++i){
    for(int j=0;j<ny;++j){
      dNy(i,j) = (y(j) == yuni(i) & ystatus(j) == 1);
      Yy(i,j) = (y(j) >= yuni(i));
    }
  }
  // Mx
  /*
  for(int i = 0; i < n1; ++i){
    for(int j=0; j < nx; ++j){
      dMx(i, j) = N(i, j) - ;
    }
  }
  */
  // Return
  return List::create(Named("dNx") = dNx, Named("dNy") = dNy, Named("lambdax") = lambdax, Named("lambday") = lambday);
}
