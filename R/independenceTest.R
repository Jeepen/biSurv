#' Non-parametric independence test for bivariate survival data
#'
#' @title Non-parametric independence test for bivariate survival data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param weight Weight function for test. Default is 'independence' which is optimal for the Frank copula 
#' @return Test statistic, SE and p-value for independence test.
#' @seealso biHazards
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
independenceTest <- function(formula, data = NULL, weight = "independence"){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
    nu1 <- nrow(eyy)
    nu2 <- ncol(eyy)
    m1 <- survfit(survival::Surv(x,xstatus)~1)
    m2 <- survfit(survival::Surv(y,ystatus)~1)
    if(class(weight) == "matrix"){}
    else if(weight == "independence"){
        weight <- m1$surv %o% m2$surv
    }
    else if(weight == "dabrowska"){
        weight <- dabrowska(formula, data)$S
    }
    else if(weight == "atRisk"){
        weight <- eyy
    }
    else if(length(weight)==1){
        weight <- matrix(weight, nrow = nu1, ncol = nu2)
    }
    else{
        stop("'weight' has to be either 'independence', 'dabrowska', 'atRisk', a number or
a matrix")
    }
    cumhaz1 <- m1$cumhaz
    cumhaz2 <- m2$cumhaz
    n <- length(x)
    M1 <- martin(x,xstatus,m1$time,cumhaz1,n,nu1)
    M2 <- martin(y,ystatus,m2$time,cumhaz2,n,nu2)
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
    out <- list(call = Call, est = testStatistic, se = sqrt(varianceEst), zval = test,
                p = 2 * pnorm(-abs(test)))
    class(out) <- "independenceTest"
    out
}

#' @export
print.independenceTest <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    coefs <- cbind(Estimate = x$est, `Std. Error` = x$se, 
                   `z value` = x$zval, `Pr(>|z|)` = x$p)
    rownames(coefs) <- "Test"
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\n")
    invisible(x)
}

#' @export
summary.independenceTest <- function(object, ...) 
{
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    coefs <- cbind(Estimate = object$est, `Std. Error` = object$se, 
                   `z value` = object$zval, `Pr(>|z|)` = object$p)
    rownames(coefs) <- "Test"
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\n")
    invisible(object)
}
