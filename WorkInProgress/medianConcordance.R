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

