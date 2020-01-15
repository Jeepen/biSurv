#' Estimator of Kendall's tau for censored data
#'
#' @title Estimator of Kendall's tau for censored data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param alpha Significance level
#' @param method Which estimator to use
#' @seealso tauPar taucpp
#' @import stats
#' @import prodlim 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
#' @useDynLib biSurv
#' @importFrom Rcpp sourceCpp
tauCens <- function(formula, data = NULL, alpha = .05, method = "adjusted"){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    id <- rep(1:length(x),2)
    if(method=="adjusted"){
        KMx <- prodlim::prodlim(prodlim::Hist(x, xstatus) ~ 1)
        KMy <- prodlim::prodlim(prodlim::Hist(y, ystatus) ~ 1)
        tauu <- taucpp(x,y,xstatus,ystatus,KMx$surv,KMy$surv,KMx$time,KMy$time)
        n <- length(x)
        a <- tauu$a
        b <- tauu$b
        var.hat <- 4 * (sum(apply(a, 2, sum)^2) - sum(a^2)) * 
            (sum(apply(b, 2, sum)^2) - sum(b^2)) / (n * (n - 1) * (n - 2)) + 
            2 * sum(a^2) * sum(b^2) / (n * (n - 1))
        var.hat <- var.hat / (sum(a^2) * sum(b^2))
        tau <- tauu$tau
    }
    else if(method == "naive"){
        obs <- (xstatus==1 & ystatus==1)
        xx <- x[obs]
        yy <- y[obs]
        tau <- cor(xx,yy,method="kendall")
        n <- length(xx)
        var.hat <- (2*(2*n+5))/(9*n*(n-1))
    }
    else{
        stop("method has to be either 'adjusted' or 'naive'")
    }
    gamma <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), data=data.frame())
    theta <- 1/exp(gamma$logtheta)
    stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame())
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    se <- c(sqrt(var.hat), sqrt(gamma$var_logtheta) * 2*exp(gamma$logtheta) / ((1 + 2*exp(gamma$logtheta))^2),
            sqrt(stable$var_logtheta) * abs(exp(stable$logtheta)/(1+exp(stable$logtheta))-exp(2*stable$logtheta)/((1+exp(stable$logtheta))^2)), NA)
    invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame())
    alpha2 <- exp(invgauss$logtheta) 
    models <- c(tauPar(theta), tauPar(alpha1,dist="posstab"), tauPar(alpha2,dist="invgauss"))
    out <- list(call = Call, coefficients = tau, se = se, models = models)
    class(out) <- "tauCens"
    out
}

#' @export
print.tauCens <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    coefficients <- c(x$coefficients, x$models)
    ## se <- c(x$se,rep(NA,3))
    se <- x$se
    coefs <- cbind(Estimate = coefficients, `Std. Error` = se, 
                   `lwr` = coefficients - 1.96 * se, `upr` = coefficients + 1.96 * se)
    rownames(coefs) <- c("Empirical", "Gamma", "Positive stable", "Inverse Gaussian")
    printCoefmat(coefs, digits = digits, na.print = "", ...)
    cat("\n")
    invisible(x)
}
