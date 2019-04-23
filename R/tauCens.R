#' Estimator of Kendall's tau for censored data
#'
#' @title Estimator of Kendall's tau for censored data
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param alpha Significance level
#' @param method Which estimator to use
#' @param without Distributions to ignore
#' @seealso tauPar taucpp
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
#' @useDynLib biSurv
#' @importFrom Rcpp sourceCpp
tauCens <- function(x,y,xstatus,ystatus, alpha = .05, method = "adjusted",without=NULL){
    models <- Names <- NULL
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
        out <- data.frame(tau = tauu$tau, SE = sqrt(var.hat),
                          lwr = tauu$tau-qnorm(1-alpha/2)*sqrt(var.hat),
                          upr = tauu$tau+qnorm(1-alpha/2)*sqrt(var.hat))
    }
    else if(method == "naive"){
        xx <- x[xstatus==1 & ystatus==1]
        yy <- y[xstatus==1 & ystatus==1]
        xxstatus <- xstatus[xstatus==1 & ystatus==1]
        yystatus <- ystatus[xstatus==1 & ystatus==1]
        tau <- cor(xx,yy,method="kendall")
        n <- length(xx)
        var.hat <- (2*(2*n+5))/(9*n*(n-1))
        out <- data.frame(tau = tau, SE = sqrt(var.hat),
                          lwr = tau-qnorm(1-alpha/2)*sqrt(var.hat),
                          upr = tau+qnorm(1-alpha/2)*sqrt(var.hat))
    }
    else{
        stop("method has to be either 'dabrowska' or 'naive'")
    }
    if(!("gamma" %in% without)){
        gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
        theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
        models <- c(models, tauPar(theta))
        Names <- c(Names, "Gamma")
    }
    if(!("posstab" %in% without)){
        stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                          control = emfrail_control(se = F, lik_ci = F, ca_test = F))
        alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
        models <- c(models, tauPar(alpha1,dist="posstab"))
        Names <- c(Names, "Positive stable")
    }
    if(!("invgauss" %in% without)){
        invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                            control = emfrail_control(se = F, lik_ci = F, ca_test = F))
        alpha2 <- exp(invgauss$logtheta)
        models <- c(models, tauPar(alpha2,dist="invgauss"))
        Names <- c(Names, "Inverse Gaussian")
    }
    list(NonParametric = out, Models = data.frame(Distribution=Names, Models=models))
}
