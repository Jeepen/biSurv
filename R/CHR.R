#' Integrated squared distance (ISD) between empirical and theoretical CHRs
#'
#' @title Integrated squared distance (ISD) between empirical and theoretical CHRs.
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data a data.frame containing the variables in the model.
#' @param n Number of steps for empirical CHR.
#' @details The CHR and the bivariate survival function are estimated seperately.
#' The ISD between the empirical CHR (as a function of the bivariate survival function)
#' and the ones implied by the different frailty models, is estimated.
#' The function returns the estimates sorted from lowest (best) to highest (worst).
#' @return Data.frame with ISD for different frailty distributions.
#' @seealso chrCpp
#' @references Chen, Min-Chi & Bandeen-Roche, Karen. (2005). A Diagnostic for Association in Bivariate Survival Models. Lifetime data analysis. 11. 245-64. 
#' @useDynLib biSurv
#' @import frailtyEM
#' @import survival
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
CHR <- function(formula, data, n = 5){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    S <- dabrowska(formula, data)
    if(length(x) != length(y)){
        stop("Length of x and y differ")
    }
    if(length(xstatus) != length(ystatus)){
        stop("Length of xstatus and ystatus differ")
    }
    n0 <- length(x)
    condis <- chrCpp(x, y, xstatus, ystatus)
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    xmin <- outer(x, x, FUN="pmin")
    ymin <- outer(y, y, FUN="pmin")
    SS <- matrix(NA, n0, n0)
    for(i in 2:n0){
        for(j in 1:(i-1)){
            SS[i,j] <- S$surv[xuni == xmin[i,j], yuni == ymin[i,j]]
        }
    }
    breaks <- seq(0,1,length.out=n+1)
    out <- numeric(n)
    for(i in 1:n){
        out[i] <- sum(condis$c[SS > breaks[i] & SS < breaks[i+1]], na.rm = T) /
            sum(condis$d[SS > breaks[i] & SS < breaks[i+1]], na.rm = T)
    }
    emp <-  rep(out, c(rep(100/n, n-1), 96-(n-1)/n*100))
    time <- c(x,y)
    status <- c(xstatus,ystatus)
    id <- rep(1:length(x),2)
    gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
    theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
    gammaCHR <- rep(1+theta, 96)
    ## gammaDiff <- mean((emp-gammaCHR)^2)
    stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                      control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    stableCHR <- 1-(1-alpha1)/(alpha1*log(seq(.01,.96,.01)))
    ## stableDiff <- mean((emp-stableCHR)^2)
    invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                        control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    invgaussCHR <- 1+1/(alpha2-log(seq(.01,.96,.01)))
    ## invgaussDiff <- mean((emp-invgaussCHR)^2)
    d <- list(call = Call, d = data.frame(S = seq(.01,.96,.01), Empirical = emp, Gamma = gammaCHR,
                                          PositiveStable = stableCHR, InverseGaussian=invgaussCHR))
    class(d) <- "CHR"
    d
}

#' @export
print.CHR <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    gammaDiff <- mean((x$d$Empirical-x$d$Gamma)^2)
    stableDiff <- mean((x$d$Empirical-x$d$PositiveStable)^2)
    invgaussDiff <- mean((x$d$Empirical-x$d$InverseGaussian)^2)
    ans <- data.frame(Distribution = c("Gamma", "Positive stable", "Inverse Gaussian"),
                      ISD = c(gammaDiff, stableDiff, invgaussDiff))
    out <- data.frame(ISD = ans$ISD[order(ans$ISD)])
    rownames(out) <- ans$Distribution[order(ans$ISD)]
    print(out)
}

#' @export
summary.CHR <- function(object, ...){
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    gammaDiff <- mean((object$d$Empirical-object$d$Gamma)^2)
    stableDiff <- mean((object$d$Empirical-object$d$PositiveStable)^2)
    invgaussDiff <- mean((object$d$Empirical-object$d$InverseGaussian)^2)
    ans <- data.frame(Distribution = c("Gamma", "Positive stable", "Inverse Gaussian"),
                      ISD = c(gammaDiff, stableDiff, invgaussDiff))
    out <- data.frame(ISD = ans$ISD[order(ans$ISD)])
    rownames(out) <- ans$Distribution[order(ans$ISD)]
    print(out)

}

#' Plot of CHR as a function of the survival function
#'
#' @title Plot of CHR as a function of the survival function.
#' @param x An object of class \code{CHR}.
#' @param ... Further arguments for \code{ggplot}.
#' @details The CHR and the bivariate survival function are estimated seperately.
#' Plotting the CHR as a function of the bivariate survival function tells us something about
#' where in the data the dependence is strongest. Is it for small failure times
#' (where the survival function is big), as implied by the positive stable model,
#' or is it for late failure times as implied by the gamma frailty model? 
#' @return Plot of CHR as a function of the survival function.
#' @seealso chrCpp
#' @references Chen, Min-Chi & Bandeen-Roche, Karen. (2005). A Diagnostic for Association in Bivariate Survival Models. Lifetime data analysis. 11. 245-64. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
plot.CHR <- function(x, ...){
    Names <- c("Gamma", "Positive stable", "Inverse Gaussian")
    melted <- melt(x$d,id="S")
    colnames(melted)[2] <- "Model"
    ggplot(melted, aes(x = .data$S, y = .data$value, color = .data$Model), ...) + geom_line() +
        geom_hline(aes(yintercept = 1), linetype = 2) + xlab(expression(S(t[1], t[2]))) +
        ylab("CHR") + theme_bw()  
}

#' Sort models according to loglikelihood
#'
#' @title Sort models according to loglikelihood.
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data a data.frame containing the variables in the model.
#' @return Data.frame with loglikelihood for different models, sorted.
#' @details It is mentioned in \code{vignette("frailtyEM_manual")} that the frailty with the
#' highest loglikelihood should be prefered since all the models are special cases of the PVF frality
#' family. This is done in a user-friendly way here. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
logLikSort <- function(formula, data){
    Call <- match.call()
    gamma <- emfrail(formula, data)
    stable <- emfrail(formula, data, distribution=emfrail_dist(dist="stable"))
    invgauss <- emfrail(formula, data, distribution=emfrail_dist(dist="pvf"))
    out <- list(call = Call, loglik = c(as.numeric(logLik(gamma)), as.numeric(logLik(stable)),
                                        as.numeric(logLik(invgauss))))
    class(out) <- "logLikSort"
    out
}

#' @export
print.logLikSort <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
                             signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    ans <- data.frame(Distribution = c("Gamma", "Positive stable", "Inverse Gaussian"),
                      loglik = x$loglik)
    out <- data.frame(loglik = rev(ans$loglik[order(ans$loglik)]))
    rownames(out) <- rev(ans$Distribution[order(ans$loglik)])
    ## printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
    ## na.print = "NA", ...)
    print(out)

}
