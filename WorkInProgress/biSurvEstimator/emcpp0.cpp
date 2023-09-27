#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix EMcpp(NumericVector x, NumericVector y, NumericVector xstatus, NumericVector ystatus, 
	  NumericVector xuni, NumericVector yuni, double tol = .001, int maxIt = 10){
  // Initialize
  double n = x.size();
  int nx = xuni.size();
  int ny = yuni.size();
  NumericMatrix p1(nx,ny);
  
  // First step
//   for(int i=0;i<nx;++i){
//     for(int j=0;j<ny;++j){
//       for(int z=0;z<n;++z){
// 	if(xstatus(z) == 1 & ystatus(z) == 1){
// 	  p1(i,j) += 1/n;	
// 	}
// 	else if(xstatus(z) == 1 & ystatus(z) == 0 & sum(yuni > y(z)) > 0){
// 	  p1(i,j) = p1(i,j) + 1/(n * sum(yuni > y(z)));
// 	}
// 	else if(xstatus(z) == 0 & ystatus(z) == 1 & sum(xuni > x(z)) > 0){
// 	  p1(i,j) = p1(i,j) + 1/(n*sum(xuni > x(z)));
// 	}
// 	else if(xstatus(z) == 0 & ystatus(z) == 0 & sum(xuni > x(z)) * sum(yuni > y(z)) > 0){
// 	  p1(i,j) = p1(i,j) + 1/(n*sum(xuni > x(z))*sum(yuni > y(z)));
// 	}
//       }
//     }
//   }
//   return p1;
// }

  for(int i=0; i<n;++i){
    if(xstatus(i) == 1 & ystatus(i) == 1){
      p1(find(xuni == x(i)), find(yuni == y(i))) = p1(find(xuni == x(i)), find(yuni == y(i))) + 
  	1/n;
    }
    else if(xstatus(i) == 1 & ystatus(i) == 0){
      if(sum(yuni > y(i)) > 0){
  	p1.elem(find(xuni == x(i)), find(yuni >= y(i))) = p1.elem(find(xuni == x(i)), find(yuni >= y(i))) + 
  	  1 / (n * sum(yuni >= y(i)));
      }
    }
    else if(xstatus(i) == 0 & ystatus(i) == 1){
      if(sum(xuni > x(i)) > 0){
  	p1.elem(find(xuni >= x(i)), find(yuni == y(i))) = p1.elem(find(xuni >= x(i)), find(yuni == y(i))) + 
  	  1 / (n * sum(xuni >= x(i)));
      }
    }    
    else{
      if(sum(xuni > x(i)) * sum(yuni > y(i)) > 0){
  	p1.elem(find(xuni >= x(i)), find(yuni >= y(i))) = p1.elem(find(xuni >= x(i)), find(yuni >= y(i))) +
  	  1 / (n * sum(xuni >= x(i)) * sum(yuni >= y(i)));
      }
    }
  }
  
  // Second step
  // NumericMatrix p2(nx,ny, fill::zeros);
  // double err = 1;
  // int itt = 0;
  // while(err > tol & itt < maxIt){
  //   itt += 1;
  //   cout << itt << endl;
  //   p2.fill(0);
  //   for(int i = 0; i < nx; ++i){
  //     for(int j = 0; j < ny; ++j){
  // 	p2(i,j) = 1/n * sum((x == xuni(i)) % (y == yuni(j)) % xstatus % ystatus);
  // 	for(int z=0;z<n;++z){
  // 	  if(sum(yuni > y(z)) > 0 & xstatus(z) == 1 & ystatus(z) == 0){
  // 	    p2(i,j) = p2(i,j) + 1/n * p1(i,j) / 
  // 	      accu(p1(find(xuni == x(z)), find(yuni > y(z))));
  // 	  }
  // 	  else if(sum(xuni > x(z)) > 0 & ystatus(z) == 1 & xstatus(z) == 0){
  // 	    p2(i,j) = p2(i,j) + 1/n * p1(i,j) / 
  // 	      accu(p1(find(xuni > x(z)), find(yuni == y(z))));	    
  // 	  }
  // 	  else if(sum(xuni > x(z))*sum(yuni > y(z)) > 0 & 
  // 		  xstatus(z) == 0 & ystatus(z) == 0){
  // 	    p2(i,j) = p2(i,j) + 1/n * p1(i,j) /
  // 	      accu(p1(find(xuni > x(z)), find(yuni > y(z))));
  // 	  }
  // 	}
  //     }
  //   }
  //   cout << p2 << endl;
  //   NumericMatrix pdiff = abs(p1-p2);
  //   err = pdiff.max();
  //   p1 = p2;
  //   // cout << p1 << endl;
  // }
// }
