#' Estimator of Kendall's tau for censored data
#'
#' @title Estimator of Kendall's tau for censored data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data a data.frame containing the variables in the model.
#' @param method which estimator to use. The non-parametric estimator from Hougaard (2000)
#' or the naive one where non-fully observed pairs are left out.
#' @seealso tauPar taucpp
#' @references Hougaard, Philip. (2000). Analysis of Multivariate Survival Data.
#' @return non-parametric estimate of Kendall's tau and parametric estimates from three different
#' frailty models.
#' @details Kendall's tau is a rank based measure of dependence. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
#' @useDynLib biSurv
#' @importFrom Rcpp sourceCpp
tauCens <- function(formula, data = NULL, method = "adjusted"){
    Call <- match.call()
    gamma <- emfrail(formula = formula, data = data)
    theta <- 1/exp(gamma$logtheta)
    stable <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="stable"))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    invgauss <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="pvf"))
    alpha2 <- exp(invgauss$logtheta) 
    models <- c(tauPar(theta), tauPar(alpha1,dist="posstab"), tauPar(alpha2,dist="invgauss"))
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    id <- rep(1:length(x),2)
    if(method=="adjusted"){
        KMx <- list(time = sort_unique(x), surv = KaplanMeier(x, xstatus))
        KMy <- list(time = sort_unique(y), surv = KaplanMeier(y, ystatus))
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
    se <- c(sqrt(var.hat), sqrt(gamma$var_logtheta) * 2*exp(gamma$logtheta) / ((1 + 2*exp(gamma$logtheta))^2),
            sqrt(stable$var_logtheta) * abs(exp(stable$logtheta)/(1+exp(stable$logtheta))-exp(2*stable$logtheta)/((1+exp(stable$logtheta))^2)),
            NA)
    out <- list(call = Call, coefficients = c(tau,models), se = se)
    class(out) <- "tauCens"
    out
}

#' @export
print.tauCens <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    coefficients <- x$coefficients
    se <- x$se
    coefs <- cbind(Estimate = coefficients, `Std. Error` = se, 
                   `lwr` = coefficients - 1.96 * se, `upr` = coefficients + 1.96 * se)
    rownames(coefs) <- c("Empirical", "Gamma", "Positive stable", "Inverse Gaussian")
    printCoefmat(coefs, digits = digits, na.print = "", ...)
    cat("\n")
    invisible(x)
}

#' Estimate of median concordance for censored data
#'
#' @title Estimate of median concordance for censored data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data a data.frame containing the variables in the model.
#' @param gamma bandwidth for pruitt estimator. Ignored if estimator is another one.
#' @param maxIt maximum number of iterations for NPMLE and pruitt estimators of
#' bivariate survival function. Doesn't do anything if estimator is 'dabrowska' (default).
#' @param method which estimator to use for the bivariate survival function.
#' @details median concordance is a rank based measure of dependence which ensures that
#' it doesn't change if we transform the marginal distributions with strictly increasing functions.
#' Median concordance is defined as
#'
#' %E(T_1,T_2)
#'
#' Median concordance is more intuitive than for instance Kendall's tau
#' since you are only comparing whether one pair of observations is concordance or
#' discordant wrt their respective medians rather than whether one pair is concordant/discordant
#' compared to another pair. Theoretically the median concordance is equal to
#' 4*S(m1,m2)-1 where m1 and m2 are the marginal medians and S is the bivariate survival function.
#' This function estimates S both non-parametrically and parametrically for three frailty models.
#' The survival function value for the frailty models at the medians are
#' uniquely determined by the copula they imply since S(m1,m2) = C(S1(m1),S2(m2)) = C(0.5,0.5),
#' where C is the survival copula implied by the model. 
#' @references Hougaard, Philip. (2000). Analysis of Multivariate Survival Data.
#' @return non-parametric estimate of median concordance and parametric estimates from
#' three different frailty models. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
medianConcordance <- function(formula, data = NULL, gamma = NULL,
                              maxIt = 100, method = "dabrowska"){
    Call <- match.call()
    gamma <- emfrail(formula = formula, data = data)
    theta <- 1/exp(gamma$logtheta)
    stable <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="stable"))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    invgauss <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="pvf"))
    alpha2 <- exp(invgauss$logtheta)
    S <- biSurv(formula, data, gamma, maxIt, method)
    gammaC <- ifelse(theta==0, 1/4, (2^(theta+1)-1)^(-1/theta))
    posstabC <- exp(-(2^alpha1*log(2)))
    invgaussC <- exp(alpha2-alpha2^(.5)*(4*(log(.5)^2/(2*alpha2)-log(.5))+alpha2)^.5)
    out <- list(Empirical = 4*S$medsurv-1, gamma = 4*gammaC-1, posstab = 4*posstabC-1,
                invgauss = 4*invgaussC-1, call = Call)
    class(out) <- "medCon"
    out
}


#' @export
print.medCon <- function (x, digits = max(3L, getOption("digits") - 3L),
                          symbolic.cor = x$symbolic.cor, ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    coefficients <- x$coefficients
    se <- x$se
    coefs <- matrix(c(x$Empirical, x$gamma, x$posstab, x$invgauss), nrow = 4, ncol = 1)
    rownames(coefs) <- c("Empirical", "Gamma", "Positive stable", 
                         "Inverse Gaussian")
    colnames(coefs) <- "Median concordance"
    printCoefmat(coefs, digits = digits, na.print = "", ...)
    cat("\n")
    invisible(x)
}


#' Get Kendall's tau from parameter or parameter from Kendall's tau
#'
#' @title Get Kendall's tau from parameter or parameter from Kendall's tau
#' @param par parameter value for frailty model.
#' @param dist frailty distribution.
#' @param output return tau or parameter.
#' @param type parameterization of frailty distribution.
#' @return Kendall's tau from parameter or parameter from Kendall's tau.
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tauPar <- function(par = 0, dist = "gamma", output = "tau", type = "alpha"){
    if(!(dist %in% c("gamma","posstab","invgauss"))){
        stop("dist has to be either 'gamma', 'posstab' or 'invgauss'")
    }
    if(!(output %in% c("tau","par"))){
        stop("output has to be either 'tau' or 'par'")
    }
    if(!(type %in% c("alpha","theta"))){
        stop("type has to be either 'alpha' or 'theta'")
    }
    if(output=="par" & (par>1|par< -1)){
        stop("'tau' has to be between -1 and 1")
    }
    if(dist=="posstab"&output=="tau"&type=="alpha"&par>=1){
        stop("alpha has to be lower than 1")
    }
    if(dist=="posstab"&output=="tau"&type=="theta"&par<=1){
        stop("theta has to be greater than 1")
    }
    if(dist=="invgauss"&output=="tau"&par<0){
        stop("'par' has to be greater than 0")
    }
    if(type=="theta"&output=="tau"&dist!="gamma"){
        par <- 1/par
    }
    if(dist == "invgauss" & par > 7e2){
        0
    }
    else{
        switch(dist, gamma={
            switch(output, tau={
                par/(2+par)
            },par={
                2*par/(1-par)
            })
        },posstab=1-par,invgauss={
            laplace <- function(alpha){
                L <- function(s) exp(alpha - sqrt(alpha) * sqrt(2*s + alpha))
                LL <- function(s){
                    exp(alpha) * (-sqrt(alpha)) * (exp(-sqrt(alpha)*sqrt(2*s+alpha))*(-sqrt(alpha)) / (2*s+alpha) -
                                                   exp(-sqrt(alpha)*sqrt(2*s+alpha))*(2*s+alpha)^(-3/2))
                }
                f <- function(s) s*L(s)*LL(s)
                4 * integrate(f, 0, Inf)$value - 1
            }
            switch(output, tau={
                laplace(par)
            },par={
                uniroot(function(alpha) laplace(alpha) - par, interval = c(.0001, 100))$root
            })
        })
    }
}

