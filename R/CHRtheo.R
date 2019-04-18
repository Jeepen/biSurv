#' Returns theoretical CHR as a function of the survival function for different frailty models
#'
#' @title Returns theoretical CHR as a function of the survival function for different frailty models
#' @param par Parameter value for frailty model
#' @param dist Frailty distribution. Gamma ("gamma"), positive stable ("posstab") and inverse Gaussian ("invgauss") are available.
#' @return A vector of length 97 with the theoretical CHR when the survival function is equal to 0.01,...,0.97
#' @seealso chrdiff CHR chr
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
CHRtheo <- function(par, dist = "gamma", type = "alpha"){
    if(!(dist %in% c("gamma","posstab","invgauss"))){
        stop("dist has to be either 'gamma', 'posstab' or 'invgauss'")
    }
    if(!(type %in% c("alpha", "theta"))){
        stop("type has to be either 'alpha' or 'theta'")
    }
    if(type=="theta" & dist != "gamma"){par<-1/par}
    switch(dist, gamma = rep(1+par, 97),
           posstab = 1 - (1-par) / (par * log(seq(.01,.97,.01))),
           invgauss=1+1/(par-log(seq(.01,.97,.01))))
}

