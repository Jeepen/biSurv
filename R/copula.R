copulaTrans <- function(t, status, km){
    indx <- match(1, km$time>t) - 1
    if(is.na(indx)) indx <- length(km$time)
    S0 <- km$surv[indx]
    if(status == 1){
        1 - S0 
    }
    else{
        time <- km$time - t
        cond <- time>0
        if(sum(cond)==0){
            1
        }
        else{
            time <- time[cond]
            survive <- km$surv[cond]
            mrl <- sum(diff(c(0,time)) * c(S0,survive[-length(survive)])) / S0
            ELT <- t + mrl
            ind <- match(1, km$time>ELT) - 1
            if(is.na(ind)) ind <- length(km$time)
            1 - km$surv[ind]
        }
    }
}

#' Plot of bivariate data, extrapolated and transformed to uniform marginals.
#'
#' @title Plot of bivariate data, extrapolated and transformed to uniform marginals.
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @return Plot of data, extrapolated, transformed to uniform marginals, and plotted
#' to compare to copulas.
#' @importFrom Rfast sort_unique
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
copulaPlot <- function(formula, data = NULL){
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    kmx <- list(time = sort_unique(x), surv = KaplanMeier(x, xstatus))
    kmy <- list(time = sort_unique(y), surv = KaplanMeier(y, ystatus))
    timex <- sapply(1:length(x), function(i){
        copulaTrans(x[i], xstatus[i], kmx)
    })
    timey <- sapply(1:length(y), function(i){
        copulaTrans(y[i], ystatus[i], kmy)
    }) 
    ## data.frame(timex,timey)
    qplot(timex, timey) + theme_bw() + xlab("") + ylab("")
}
