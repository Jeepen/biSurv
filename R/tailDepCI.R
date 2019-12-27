#' Bootstrapped confidence interval for "tail dependence"
#'
#' @title Bootstrapped confidence interval for "tail dependence"
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param cluster Cluster variable
#' @param q Quantile to estimate "tail-dependence" for
#' @param method How to estimate survival probabilities. Available are 'dabrowska' and 'fast'
#' @param tail Tail to estimate "tail dependence" for
#' @param n Number of bootstraps
#' @return CI for "tail dependence"
#' @seealso tailDep
#' @import boot
#' @import graphics
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDepCI <- function(formula, data, cluster, q, method = "dabrowska", tail="lwr", n = 1000){
    d <- uniTrans(formula, data, cluster)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    d <- data.frame(x=x,y=y,xstatus=xstatus,ystatus=ystatus)
    f <- function(data,i){
        d <- data[i,]
        tailDependence(d$x,d$y,d$xstatus,d$ystatus,q,method=method,tail=tail,without=c("gamma","posstab","invgauss"))[1,2]
    }
    out <- boot(data=d,statistic=f,R=n)
    print(boot.ci(out))
}
