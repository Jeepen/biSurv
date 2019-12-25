#' Returns non-parametric estimates of the bivariate hazard function
#'
#' @title Returns non-parametric estimates of the bivariate hazard function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame in which to interpret the variables named in the
#'          \code{formula} argument.
#' @return List with the estimates of the three bivariate hazard functions. Each hazard function is a matrix.
#' @seealso hazardscpp dabrowska
#' @import Rfast
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
biHazards <- function(formula, data){
    d <- uniTrans(formula, data)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    h <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
    h$lambda11[is.nan(h$lambda11)] <- 0
    h$lambda10[is.nan(h$lambda10)] <- 0
    h$lambda01[is.nan(h$lambda01)] <- 0
    h
}
