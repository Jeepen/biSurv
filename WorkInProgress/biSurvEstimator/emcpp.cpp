#include "RcppArmadillo.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
mat EMcpp(const vec &x, const vec &y, const vec &xstatus, const vec &ystatus, 
	  const vec &xuni, const vec &yuni, double tol = .001, int maxIt = 10){
  // Initialize
  double n = x.size();
  int nx = xuni.size();
  int ny = yuni.size();
  mat p1(nx,ny, fill::zeros);
  // mat p3(nx,ny, fill::zeros);
  
  // First step
  for(int i=0; i<n;++i){
    if(xstatus(i) == 1 & ystatus(i) == 1){
      p1(find(xuni == x(i)), find(yuni == y(i))) += 1/n;
    }
    else if(xstatus(i) == 1 & ystatus(i) == 0){
      if(sum(yuni > y(i)) > 0){
	p1(find(xuni == x(i)), find(yuni > y(i))) += 1/(n * sum(yuni > y(i)));
      }
    }
    else if(xstatus(i) == 0 & ystatus(i) == 1){
      if(sum(xuni > x(i)) > 0){
	p1(find(xuni > x(i)), find(yuni == y(i))) += 1/(n * sum(xuni > x(i)));
      }
    }    
    else{
      if(sum(xuni > x(i)) * sum(yuni > y(i)) > 0){
	p1(find(xuni > x(i)), find(yuni > y(i))) += 1/(n * sum(xuni > x(i)) * sum(yuni > y(i)));
      }
    }
  }
  // p1 = p1 + p3;
  
  // Second step
  // mat p2(nx,ny, fill::zeros);
  // mat p2 = p3;
  mat p2(nx,ny, fill::zeros);
  mat pdiff(nx,ny);
  double err = 1;
  int itt = 0;
  while(err > tol & itt < maxIt){
    itt += 1;
    cout << itt << endl;
    // p2 = p3;
    p2.zeros();
    // cout << p2.max() << endl;
    for(int i = 0; i < nx; ++i){
      for(int j = 0; j < ny; ++j){
  	p2(i,j) = 1/n * sum((x == xuni(i)) % (y == yuni(j)) % xstatus % ystatus);
  	for(int z=0;z<n;++z){
  	  if(xstatus(z) == 1 & ystatus(z) == 0 & xuni(i) == x(z) & yuni(j) > y(z)){
  	    p2(i,j) += 1/n * p1(i,j) / accu(p1(find(xuni == x(z)), find(yuni > y(z))));
  	  }
  	  else if(ystatus(z) == 1 & xstatus(z) == 0 & xuni(i) > x(z) & yuni(j) == y(z)){
  	    p2(i,j) += 1/n * p1(i,j) / accu(p1(find(xuni > x(z)), find(yuni == y(z))));	    
  	  }
  	  else if(xstatus(z) == 0 & ystatus(z) == 0 & xuni(i) > x(z) & yuni(j) > y(z)){
  	    p2(i,j) += 1/n * p1(i,j) / accu(p1(find(xuni > x(z)), find(yuni > y(z))));
  	  }
  	}
      }
    }
    // cout << p2 << endl;
    pdiff = abs(p1-p2);
    err = pdiff.max();
    cout << err << endl;
    p1 = p2;
    // cout << p1 << endl;
  }
  return p1;
}
