#' Returns non-parametric estimates of the bivariate hazard function
#'
#' @title Returns non-parametric estimates of the bivariate hazard function
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @return List with the estimates of the three bivariate hazard functions. Each hazard function is a matrix.
#' @seealso hazardscpp dabrowska
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
biHazards <- function(x,y,xstatus,ystatus){
  xuni <- sort_unique(x)
  yuni <- sort_unique(y)
  h <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
  h$lambda11[is.nan(h$lambda11)] <- 0
  h$lambda10[is.nan(h$lambda10)] <- 0
  h$lambda01[is.nan(h$lambda01)] <- 0
  h
}
