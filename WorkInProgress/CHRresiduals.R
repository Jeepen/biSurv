CHRres <- function(x,y,xstatus,ystatus, weight = "auto"){
    h <- biHazards(x,y,xstatus,ystatus)
    if(weight == "auto"){
        weight <- eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
    }
    else if(all(dim(weight) != dim(h$lambda11))){
        stop("dim(weight) has to equal dim(lambda11)")
    }
    else{
        eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
    }
    CHRresiduals <- weight * (h$lambda11 - h$lambda10 * h$lambda01)
    CHRresiduals[is.nan(CHRresiduals)] <- 0
    testStatistic <- sum(CHRresiduals)
    m1 <- survfit(Surv(x,xstatus)~1)
    m2 <- survfit(Surv(y,ystatus)~1)
    ## Variance estimate
    help <- c(m1$cumhaz[1],diff(m1$cumhaz)) %o% c(m2$cumhaz[1],diff(m2$cumhaz))
    varianceEst <- sum(weight^2 * eyy * help)
    ## varianceEst <- 
    ## Present wrapping
    test <- -abs(testStatistic / sqrt(varianceEst))
    d <- data.frame(testStatistic, sqrt(varianceEst), 2*pnorm(test))
    colnames(d) <- c("Test statistic", "SE", "p-value")
    d
}
