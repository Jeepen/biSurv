#' Simulate data from frailty model with Weibull baseline hazard
#'
#' @title Simulate data from frailty model with Weibull baseline hazard
#' @param n Number of observations in each cluster.
#' @param dist Frailty distribution.
#' @param par Frailty parameter.
#' @param shape shape parameter of Weibull distribution.
#' @param scale scale parameter of Weibull distribution.
#' @param cov1 Covariates for the first vector
#' @param cov2 Covariates for the second vector
#' @param beta Conditional (log-) hazard ratios
#' @importFrom stabledist rstable
#' @importFrom statmod rinvgauss
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
simFrail <- function(n = 100, dist = "gamma", par = 1, shape = 1, scale = 1, cov1 = 1, cov2 = 1, beta = 0){
    u1 <- runif(n)
    u2 <- runif(n)
    if(dist %in% c("gamma", "invgauss") & par == 0) U <- 1
    else if(dist == "gamma") U <- rgamma(n, 1/par, 1/par)
    else if(dist == "posstab") U <- stabledist::rstable(n, par, beta = 1, pm = 1)
    else if(dist == "invgauss") U <- statmod::rinvgauss(n, shape = par)
    if(is.null(dim(cov1))){
        x <- (-log(1-u1) / (scale * exp(beta * cov1) * U))^(1/shape)
    }
    if(is.null(dim(cov2))){
        y <- (-log(1-u2) / (scale * exp(beta * cov2) * U))^(1/shape)
    }
    if(length(cov1) == n) data.frame(id = rep(1:n,2), time = c(x,y), status = 1, cov = c(cov1,cov2))
    else data.frame(id = rep(1:n,2), time = c(x,y), status = 1)
}






