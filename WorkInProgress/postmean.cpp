#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List NY(NumericVector time, NumericVector status, NumericVector timeuni, int nuni, int n){
  NumericMatrix N(nuni,n);
  NumericMatrix Y(nuni,n);
  for(int t = 0; t<nuni; ++t){
    for(int i = 0; i<n;++i){
      Y(t,i) = time(i) >= timeuni(t);
      N(t,i) = status(i) * (time(i) <= timeuni(t));
    }
  }
  return List::create(Named("N") = N, Named("Y") = Y);
}

// [[Rcpp::export]]
cube gghelper(NumericMatrix Lambda, NumericVector time, 
              NumericVector LP, NumericMatrix R, NumericMatrix xi, 
              NumericMatrix Y, double theta, int n0, int n1){
  cube g = cube(n0, n0, n1);
  for(int u = 0; u < n0; ++u){
    std::cout << u << std::endl;
    for(int t = 0; t < n0;++t){
      for(int i = 0; i < n1; ++i){
        g(u, t, i) = theta / R(t,i) * xi(t,i) * 
          (exp(theta * Lambda(t,2*i) * exp(LP(2*i))) *
          exp(LP(2*i)) * Y(u,2*i) +
          exp(theta * Lambda(t,2*i+1) * exp(LP(2*i+1))) *
          exp(LP(2*i+1)) * Y(u,2*i+1));
      }
    }
  }
  return g;
}

// [[Rcpp::export]]
NumericMatrix pihelper(NumericVector LP, NumericVector H, NumericVector RR, double theta,
		       NumericMatrix Ndot, NumericMatrix Y, NumericMatrix N, int nuni, int n){
  NumericMatrix out(nuni,n/2);
  for(int t=0; t<nuni; ++t){
    for(int i=0; i<(n/2); ++i){
      out(t,i) = exp(LP(2*i)) * Y(t,2*i) * 
	(-(1/theta + Ndot(nuni,i)) / RR(i) * (1 + theta * H(2*i)) *
	 exp(theta * H(2*i)) + exp(theta * H(2*i)) / (theta * RR(i)) +
	 N(nuni,2*i) + (1 + theta * Ndot[nuni,i]) * exp(theta * H(2*i)) /
	 pow(RR(i),2) * (H(2*i) * exp(theta * H(2*i)) + 
		    H(2*i+1) * exp(theta * H(2*i+1)))) +
	exp(LP(2*i+1)) * Y(t,2*i+1) * 
	(-(1/theta + Ndot(nuni,i)) / RR(i) * (1 + theta * H(2*i+1)) *
	 exp(theta * H(2*i+1)) + exp(theta * H(2*i+1)) / (theta * RR(i)) +
	 N(nuni,2*i+1) + (1 + theta * Ndot[nuni,i]) * exp(theta * H(2*i+1)) /
	 pow(RR(i),2) * (H(2*i) * exp(theta * H(2*i)) + 
		    H(2*i+1) * exp(theta * H(2*i+1))));
    }
  }
  return out;
}

// [[Rcpp::export]]
List SEVcpp(NumericVector lp, NumericMatrix Y, NumericMatrix X, int nuni, int n, int colX){
  NumericVector S0(nuni);
  NumericMatrix S1(nuni, colX);
  cube S2(nuni, colX, colX);
  for(int t=0; t<nuni;++t){
    S0(t) = sum(Y(t,_) * exp(lp)) / (n/2);
    S1(t,_) = sum(Y(t,_) * exp(lp) * X(i,_)
    for(int i=0; i<colX; ++i){
      S1(t,i) = sum(Y(t,_) * exp(lp) * X(i,_))/ (n/2);
    }
  }
  return List::create(Named("S0") = S0, Named("S1") = S1); 
}
