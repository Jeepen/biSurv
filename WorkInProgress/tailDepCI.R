    
#' Bootstrapped confidence interval for "tail dependence"
#'
#' @title Bootstrapped confidence interval for "tail dependence"
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param q Quantile to estimate "tail-dependence" for
#' @param method How to estimate survival probabilities. Available are 'dabrowska' and 'fast'
#' @param tail Tail to estimate "tail dependence" for
#' @param n Number of bootstraps
#' @param level The confidence level required.
#' @param type Type of confidence interval. Available are "normal", which is equal to estimate plus/minus 1.96 times the bootstrap standard error,
#' and "quantile", which is the empirical quantiles of the bootstrapped estimates. 
#' @return CI for "tail dependence"
#' @seealso tailDep
#' @import graphics
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDepCI <- function(formula, data, q, method = "fast", tail="lwr", n = 1000, level = .95, type = "normal"){
    if(!(type %in% c("normal", "quantile")))
        stop("'type' has to be either 'normal' or 'quantile'")
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    data <- data.frame(time = c(x,y), status = c(xstatus,ystatus), cluster = rep(1:length(x),2))
    sims <- numeric(n)
    for(i in 1:n){
        ind <- sample(1:length(x), replace = TRUE)
        d <- data.frame(time = c(x[ind], y[ind]), status = c(xstatus[ind], ystatus[ind]),
                        id = rep(1:length(x),2))
        sims[i] <- tailDep(Surv(time,status) ~ cluster(id), data = d, q, method = method, tail=tail)
    }
    d <- data.frame(time = c(x, y), status = c(xstatus, ystatus), id = rep(1:length(x),2))
    tmp <- tailDep(Surv(time,status) ~ cluster(id), data = d, q, method = method, tail=tail)
    a <- (1-level)/2
    a <- c(a,1-a)
    ## pct <- format.perc(a, 3)
    pct <- paste(a*100, "%")
    ci <- array(NA_real_, dim = c(1L, 2L), dimnames = list("Tail dependence", pct))  
    if(type=="normal"){
        ci[] <- tmp +  qnorm(a) * sd(sims)
    }
    else{
        ci[] <- quantile(a, sims)
    }
    ci
}
