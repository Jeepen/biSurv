#include <Rcpp.h>
using namespace Rcpp;

//' Concordance and discordance matrices for bivariate survival data
//' 
//' @title Concordance and discordance matrices for bivariate survival data
//' @param x,y Vectors of failure times
//' @param xstatus,ystatus Status indicators for failure times
//' @return Concordance and discordance matrices for bivariate survival data
//' @seealso CHR CHRtheo chrdiff
//' @export
//' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
// [[Rcpp::export]]
List chrCpp(NumericVector x, NumericVector y, NumericVector xstatus, NumericVector ystatus){
  double n = x.size();
  NumericMatrix c = NumericMatrix(n);
  NumericMatrix d = NumericMatrix(n);
  for(int i=0; i<n; ++i){
    for(int j=0; j<n; ++j){
      if((xstatus(i) == 1) & (ystatus(i) == 1) & 
         (xstatus(j) == 1) & (ystatus(j) == 1)){
        c(i,j) = ((x(i) > x(j)) & (y(i) > y(j))) + 
          ((x(i) < x(j)) && (y(i) < y(j)));
        d(i,j) = ((x(i) > x(j)) & (y(i) < y(j))) + 
          ((x(i) < x(j)) & (y(i) > y(j)));
      }
      else if((xstatus(i) == 0) & (ystatus(i) == 0) &
              (xstatus(j) == 1) & (ystatus(j) == 1)){
        c(i,j) = ((x(i) > x(j)) & (y(i) > y(j)));
      }
      else if((xstatus(i) == 1) & (ystatus(i) == 1) &
              (xstatus(j) == 0) & (ystatus(j) == 0)){
        c(i,j) = ((x(i) < x(j)) & (y(i) < y(j)));
      }
      else if((xstatus(i) == 1) & (ystatus(i) == 0) &
              (xstatus(j) == 0) & (ystatus(j) == 1)){
        d(i,j) = ((x(i) < x(j)) & (y(i) > y(j)));
      }
      else if((xstatus(i) == 0) & (ystatus(i) == 1) &
              (xstatus(j) == 1) & (ystatus(j) == 0)){
        d(i,j) = ((x(i) > x(j)) & (y(i) < y(j)));
      }
      else if((xstatus(i) == 1) & (ystatus(i) == 1) &
              (xstatus(j) == 1) & (ystatus(j) == 0)){
        c(i,j) = ((x(i) < x(j)) & (y(i) < y(j)));
        d(i,j) = ((x(i) > x(j)) & (y(i) < y(j)));
      }
      else if((xstatus(i) == 1) & (ystatus(i) == 0) &
              (xstatus(j) == 1) & (ystatus(j) == 1)){
        c(i,j) = ((x(i) > x(j)) & (y(i) > y(j)));
        d(i,j) = ((x(i) < x(j)) & (y(i) > y(j)));
      }
      else if((xstatus(i) == 0) & (ystatus(i) == 1) &
              (xstatus(j) == 1) & (ystatus(j) == 1)){
        c(i,j) = ((x(i) > x(j)) & (y(i) > y(j)));
        d(i,j) = ((x(i) > x(j)) & (y(i) < y(j)));
      }
      else if((xstatus(i) == 1) & (ystatus(i) == 1) &
              (xstatus(j) == 0) & (ystatus(j) == 1)){
        c(i,j) = ((x(i) < x(j)) & (y(i) < y(j)));
        d(i,j) = ((x(i) < x(j)) & (y(i) > y(j)));
      }
    }
  }
  return List::create(Named("c") = c, Named("d") = d);
}
