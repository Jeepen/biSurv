#' Returns non-parametric estimate of piecewise constant CHR
#'
#' @title Returns non-parametric estimate of piecewise constant CHR
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param n Number of 'pieces' to estimate CHR on
#' @return Non-parametric estimate of piecewise constant CHR
#' @seealso CHRtheo chrdiff chrCpp
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
CHR <- function(x,y,xstatus,ystatus,n=5){
  n0 <- length(x)
  S <- dabrowska(x,y,xstatus,ystatus)
  condis <- chrCpp(x, y, xstatus, ystatus)
  xuni <- sort_unique(x)
  yuni <- sort_unique(y)
  xmin <- outer(x, x, FUN = "pmin")
  ymin <- outer(y, y, FUN = "pmin")
  SS <- matrix(NA,n0,n0)
  for(i in 2:n0){
    for(j in 1:(i-1)){
      SS[i,j] <- S$S[xuni == xmin[i,j], yuni == ymin[i,j]]
    }
  }
  breaks <- seq(0,1,length.out=n+1)
  out <- numeric(n)
  for(i in 1:n){
      out[i] <- sum(condis$c[SS > breaks[i] & SS < breaks[i+1]], na.rm = T) /
          sum(condis$d[SS > breaks[i] & SS < breaks[i+1]], na.rm = T)
  }
  out
}
