#' Get Kendall's tau from parameter or parameter from Kendall's tau
#'
#' @title Get Kendall's tau from parameter or parameter from Kendall's tau
#' @param par Parameter for frailty model
#' @param dist Frailty distribution
#' @param output Return tau or parameter
#' @param type Parameterization of frailty distribution
#' @return Kendall's tau from parameter or parameter from Kendall's tau
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tauPar <- function(par = 0, dist = "gamma", output = "tau", type = "alpha"){
    if(!(dist %in% c("gamma","posstab","invgauss"))){
        stop("dist has to be either 'gamma', 'posstab' or 'invgauss'")
    }
    if(!(output %in% c("tau","par"))){
        stop("output has to be either 'tau' or 'par'")
    }
    if(!(type %in% c("alpha","theta"))){
        stop("type has to be either 'alpha' or 'theta'")
    }
    if(type=="theta"&output=="tau"&dist!="gamma"){
        par <- 1/par
    }
    switch(dist, gamma={
        switch(output, tau={
            par/(2+par)
        },par={
            2*par/(1-par)
        })
    },posstab=1-par,invgauss={
        laplace <- function(alpha){
            L <- function(s) exp(alpha - sqrt(alpha) * sqrt(2*s + alpha))
            LL <- function(s){
                exp(alpha) * (-sqrt(alpha)) * (exp(-sqrt(alpha)*sqrt(2*s+alpha))*(-sqrt(alpha)) / (2*s+alpha) -
                                               exp(-sqrt(alpha)*sqrt(2*s+alpha))*(2*s+alpha)^(-3/2))
            }
            f <- function(s) s*L(s)*LL(s)
            4 * integrate(f, 0, Inf)$value - 1
        }
        switch(output, tau={
            laplace(par)
        },par={
            uniroot(function(alpha) laplace(alpha) - par, interval = c(.0001, 100))$root
        })
    })
}
