#' Estimator of Kendall's tau for censored data
#'
#' @title Estimator of Kendall's tau for censored data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame in which to interpret the variables named in the
#'          \code{formula} argument.
#' @param alpha Significance level
#' @param method Which estimator to use
#' @seealso tauPar taucpp
#' @import stats
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
#' @useDynLib biSurv
#' @importFrom Rcpp sourceCpp
tauCens <- function(formula, data, alpha = .05, method = "adjusted"){
    Call <- match.call()
    models <- NULL
    d <- uniTrans(formula, data)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    id <- rep(1:length(x),2)
    if(method=="adjusted"){
        KMx <- prodlim::prodlim(Hist(x, xstatus) ~ 1)
        KMy <- prodlim::prodlim(Hist(y, ystatus) ~ 1)
        tauu <- taucpp(x,y,xstatus,ystatus,KMx$surv,KMy$surv,KMx$time,KMy$time)
        n <- length(x)
        a <- tauu$a
        b <- tauu$b
        var.hat <- 4 * (sum(apply(a, 2, sum)^2) - sum(a^2)) * 
            (sum(apply(b, 2, sum)^2) - sum(b^2)) / (n * (n - 1) * (n - 2)) + 
            2 * sum(a^2) * sum(b^2) / (n * (n - 1))
        var.hat <- var.hat / (sum(a^2) * sum(b^2))
        tau <- tauu$tau
        ## out <- data.frame(tau = tauu$tau, SE = sqrt(var.hat),
        ##                   lwr = tauu$tau-qnorm(1-alpha/2)*sqrt(var.hat),
        ##                   upr = tauu$tau+qnorm(1-alpha/2)*sqrt(var.hat))
    }
    else if(method == "naive"){
        xx <- x[xstatus==1 & ystatus==1]
        yy <- y[xstatus==1 & ystatus==1]
        xxstatus <- xstatus[xstatus==1 & ystatus==1]
        yystatus <- ystatus[xstatus==1 & ystatus==1]
        tau <- cor(xx,yy,method="kendall")
        n <- length(xx)
        var.hat <- (2*(2*n+5))/(9*n*(n-1))
    }
    else{
        stop("method has to be either 'adjusted' or 'naive'")
    }
    gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
    theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1] 
    stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                      control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                        control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    models <- c(tauPar(theta), tauPar(alpha1,dist="posstab"), tauPar(alpha2,dist="invgauss"))
    list(call = Call, coefficients = tau, se = sqrt(var.hat), models = models)
}

#' @export
print.tauCens <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    2
}
