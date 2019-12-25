#' Non-parametric independence test for bivariate survival data
#'
#' @title Non-parametric independence test for bivariate survival data
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param weight Weight function for test. Default is 'independence' which is optimal for the Frank copula
#' @return Test statistic, SE and p-value for independence test.
#' @seealso biHazards
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
independenceTest <- function(x,y,xstatus,ystatus, weight = "independence"){
    Call <- match.call()
    eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
    nu1 <- nrow(eyy)
    nu2 <- ncol(eyy)
    m1 <- survfit(Surv(x,xstatus)~1)
    m2 <- survfit(Surv(y,ystatus)~1)
    if(class(weight) == "matrix"){}
    else if(weight == "independence"){
        weight <- m1$surv %o% m2$surv
    }
    else if(weight == "dabrowska"){
        weight <- dabrowska(x,y,xstatus,ystatus)$S
    }
    else if(weight == "atRisk"){
        weight <- eyy
    }
    else if(length(weight)==1){
        weight <- matrix(weight, nrow = nu1, ncol = nu2)
    }
    else{
        stop("'weight' has to be either 'independence', 'dabrowska', 'atRisk', a number or
a specific matrix")
    }
    cumhaz1 <- m1$cumhaz
    cumhaz2 <- m2$cumhaz
    n <- length(x)
    M1 <- matrix(NA, nrow = nu1, ncol = n)
    M2 <- matrix(NA, nrow = nu2, ncol = n)
    for(t in 1:nu1){
        for(i in 1:n){
            M1[t,i] <- xstatus[i] * (x[i] <= m1$time[t]) - cumhaz1[t]*(x[i] > m1$time[t]) -
                cumhaz1[which(m1$time == x[i])]*(x[i] <= m1$time[t])
        }
    }
    for(t in 1:nu2){
        for(i in 1:n){
            M2[t,i] <- ystatus[i] * (y[i] <= m2$time[t]) - cumhaz2[t]*(y[i] > m2$time[t]) -
                cumhaz2[which(m2$time == y[i])]*(y[i] <= m2$time[t])
        }
    }
    dM1 <- rbind(M1[1,],diff(M1))
    dM2 <- rbind(M2[1,],diff(M2))
    testStatistic <- 0
    for(i in 1:n){
       testStatistic <- testStatistic + sum((dM1[,i] %o% dM2[,i]) * weight)
    }
    ## Variance estimate
    help <- c(cumhaz1[1],diff(cumhaz1)) %o% c(cumhaz2[1],diff(cumhaz2))
    varianceEst <- sum(weight^2 * eyy * help) * n
    ## Present wrapping
    test <- testStatistic / sqrt(varianceEst)
    ## d <- data.frame(testStatistic, sqrt(varianceEst), test, 2*pnorm(-abs(test)))
    ## colnames(d) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    out <- list(call = Call, est = testStatistic, se = sqrt(varianceEst), zval = test,
                p = 2 * pnorm(-abs(test)))
    ## out$call <- Call
    ## out$coefficients <- 
    class(out) <- "independenceTest"
    out
}

#' @export
print.independenceTest <- function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    ## else cat("\nCoefficients:\n")
    coefs <- cbind(Estimate = x$est, `Std. Error` = x$se, 
                   `z value` = x$zval, `Pr(>|t|)` = x$p)
    rownames(coefs) <- "Test"
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\n")
    invisible(x)
}
