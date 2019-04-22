#' Bootstrapped confidence interval for "tail dependence"
#'
#' @title Bootstrapped confidence interval for "tail dependence"
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param q Quantile to estimate "tail-dependence" for
#' @param method How to estimate survival probabilities. Available are 'dabrowska' and 'fast'
#' @param tail Tail to estimate "tail dependence" for
#' @param n Number of bootstraps
#' @return CI for "tail dependence"
#' @seealso tailDep
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDepCI <- function(x,y,xstatus,ystatus,q,method="dabrowska",tail="lwr",n=1000){
    d <- data.frame(x=x,y=y,xstatus=xstatus,ystatus=ystatus)
    f <- function(data,i){
        d <- data[i,]
        tailDep(d$x,d$y,d$xstatus,d$ystatus,q,method=method,tail=tail,without=c("gamma","posstab","invgauss"))[1,2]
    }
    out <- boot::boot(data=d,statistic=f,R=n)
    print(boot::boot.ci(out))
    plot(out)
}
