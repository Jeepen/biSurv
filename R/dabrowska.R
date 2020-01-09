#' Dabrowska estimator of the bivariate survival function
#'
#' @title Dabrowska estimator of the bivariate survival function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @return Matrix with Dabrowska estimate of the bivariate survival function
#' @seealso biHazards hazardscpp
#' @import prodlim
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
dabrowska <- function(formula, data){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
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
