#' Integrated squared distance (ISD) between empirical and theoretical CHRs
#'
#' @title Integrated squared distance (ISD) between empirical and theoretical CHRs
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param n Number of steps for empirical CHR
#' @return Data.frame with ISD for different frailty distributions
#' @seealso chrCpp
#' @references Chen, Min-Chi & Bandeen-Roche, Karen. (2005). A Diagnostic for Association in Bivariate Survival Models. Lifetime data analysis. 11. 245-64. 
#' @useDynLib biSurv
#' @import frailtyEM
#' @import survival
#' @import ggplot2
#' @import reshape2
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
CHR <- function(formula, data, n = 5){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    S <- dabrowska(formula, data)
    if(length(x)!=length(y)){
        stop("Length of x and y differ")
    }
    if(length(xstatus)!=length(ystatus)){
        stop("Length of xstatus and ystatus differ")
    }
    n0 <- length(x)
    condis <- chrCpp(x,y,xstatus,ystatus)
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    xmin <- outer(x,x,FUN="pmin")
    ymin <- outer(y,y,FUN="pmin")
    SS <- matrix(NA,n0,n0)
    for(i in 2:n0){
        for(j in 1:(i-1)){
            SS[i,j] <- S$S[xuni == xmin[i,j], yuni == ymin[i,j]]
        }
    }
    breaks <- seq(0,1,length.out=n+1)
    out <- numeric(n)
    for(i in 1:n){
        out[i] <- sum(condis$c[SS > breaks[i] & SS < breaks[i+1]], na.rm = T) /
            sum(condis$d[SS > breaks[i] & SS < breaks[i+1]], na.rm = T)
    }
    emp <-  rep(out, c(rep(100/n, n-1), 97-(n-1)/n*100))
    time <- c(x,y)
    status <- c(xstatus,ystatus)
    id <- rep(1:length(x),2)
    gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
    theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
    gammaCHR <- rep(1+theta, 97)
    ## gammaDiff <- mean((emp-gammaCHR)^2)
    stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                      control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    stableCHR <- 1-(1-alpha1)/(alpha1*log(seq(.01,.97,.01)))
    ## stableDiff <- mean((emp-stableCHR)^2)
    invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                        control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    invgaussCHR <- 1+1/(alpha2-log(seq(.01,.97,.01)))
    ## invgaussDiff <- mean((emp-invgaussCHR)^2)
    d <- list(call = Call, d = data.frame(S = seq(.01,.97,.01), Empirical = emp, Gamma = gammaCHR,
                                          PositiveStable = stableCHR, InverseGaussian=invgaussCHR))
    class(d) <- "CHR"
    d
}

#' @export
print.CHR <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
    signif.stars = getOption("show.signif.stars"), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    gammaDiff <- mean((x$d$Empirical-x$d$Gamma)^2)
    stableDiff <- mean((x$d$Empirical-x$d$PositiveStable)^2)
    invgaussDiff <- mean((x$d$Empirical-x$d$InverseGaussian)^2)
    ans <- data.frame(Distribution = c("Gamma", "Positive stable", "Inverse Gaussian"),
                      ISD = c(gammaDiff, stableDiff, invgaussDiff))
    out <- data.frame(ISD = ans$ISD[order(ans$ISD)])
    rownames(out) <- ans$Distribution[order(ans$ISD)]
    ## printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
    ## na.print = "NA", ...)
    print(out)
    ## cat("\n")
    ## invisible(x)
}

#' @export
summary.CHR <- function(object, ...) 
{
    cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    gammaDiff <- mean((object$d$Empirical-object$d$Gamma)^2)
    stableDiff <- mean((object$d$Empirical-object$d$PositiveStable)^2)
    invgaussDiff <- mean((object$d$Empirical-object$d$InverseGaussian)^2)
    ans <- data.frame(Distribution = c("Gamma", "Positive stable", "Inverse Gaussian"),
                      ISD = c(gammaDiff, stableDiff, invgaussDiff))
    out <- data.frame(ISD = ans$ISD[order(ans$ISD)])
    rownames(out) <- ans$Distribution[order(ans$ISD)]
    ## printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
    ## na.print = "NA", ...)
    print(out)
    ## cat("\n")
    ## invisible(x)
}



#' Plot of CHR as a function of the survival function
#'
#' @title Plot of CHR as a function of the survival function
#' @param x An object of class \code{CHR}
#' @param ... Further arguments for \code{ggplot}
#' @return Plot of CHR as a function of the survival function
#' @seealso chrCpp
#' @references Chen, Min-Chi & Bandeen-Roche, Karen. (2005). A Diagnostic for Association in Bivariate Survival Models. Lifetime data analysis. 11. 245-64. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
plot.CHR <- function(x, ...){
    Names <- c("Gamma", "Positive stable", "Inverse Gaussian")
    melted <- melt(x$d,id="S")
    colnames(melted)[2] <- "Model"
    ggplot(melted, aes(x=S,y=value,color=Model),...) + geom_line() +
        geom_hline(aes(yintercept=1),linetype=2) + xlab(expression(S(t[1],t[2]))) +
        ylab("CHR") + theme_bw()  
}
