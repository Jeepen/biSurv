#' Get integrated squared difference (ISD) between empirical and theoretical CHR
#'
#' @title Get integrated squared difference (ISD) between empirical and theoretical CHR
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param par Parameter value for frailty model
#' @param dist Frailty distribution. Gamma ("gamma"), positive stable ("posstab") and inverse Gaussian ("invgauss") are available.
#' @param type Parameterization of frailty parameter (either "alpha" or "theta")
#' @param n How many intervals should the CHR be estimated on.
#' @return Estimated ISD between theoretical and empirical CHR
#' @seealso CHRtheo CHR chrCpp
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
chrdiff <- function(x, y, xstatus, ystatus, par==NULL, dist=NULL, type = "alpha", n = 5){
  imp <- CHR(x,y,xstatus,ystatus,n=n)
  imp0 <- numeric(97)
  breaks <- seq(0,1,length.out=n+1)
  for(i in 1:(n-1)){
    imp0[(1+(i-1)*100/n):(i*100/n)] <- imp[i]
  }
  imp0[(1+(n-1)*100/n):97] <- imp[n]
  theo <- CHRtheo(par, dist = dist, type = type)
  mean((imp0 - theo)^2)
}
