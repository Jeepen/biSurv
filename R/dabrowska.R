#' Dabrowska estimator of the bivariate survival function
#'
#' @title Dabrowska estimator of the bivariate survival function
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @return Matrix with Dabrowska estimate of the bivariate survival function
#' @seealso hazards hazardscpp
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
dabrowska <- function(x,y,xstatus,ystatus){
  xuni <- sort_unique(x)
  yuni <- sort_unique(y)
  haz <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
  KMx <- prodlim(Hist(x,xstatus) ~ 1)
  KMy <- prodlim(Hist(y,ystatus) ~ 1)
  H <- (haz$lambda10*haz$lambda01-haz$lambda11) / 
    ((1-haz$lambda10)*(1-haz$lambda01))
  H[is.nan(H)] <- 0
  list(S = t(apply(apply(1 - H, 2, cumprod), 1, cumprod)) * 
         (KMx$surv %*% t(KMy$surv)), timex = xuni, timey = yuni)
}
