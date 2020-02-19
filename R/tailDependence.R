#' Estimators of "tail dependence" for bivariate survival data
#'
#' @title Estimators of "tail-dependence" for bivariate survival data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param q Quantile to estimate "tail-dependence" for
#' @param tail Tail to estimate "tail dependence" for
#' @param method What estimator to use
#' @return Estimate of "tail-dependence"
#' @seealso tailDepCI tailDependencePlot
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDependence <- function(formula, data, q, tail = "lwr", method = "fast"){
    Call <- match.call()
    if(!(method %in% c("dabrowska","fast"))){
        stop("method has to be either 'dabrowska' or 'fast'")
    }
    if(!(tail %in% c("lwr","upr"))){
        stop("tail has to be either 'lwr' or 'upr'")
    }
    ## Models
    gamma <- emfrail(formula = formula, data = data, control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    ## theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
    theta <- 1 / exp(gamma$logtheta)
    out <- switch(tail,lwr=((2*(1-q)^(-theta) - 1)^(-1/theta) + 2*q - 1) / q,
                        upr=(2*(1-q)^(-theta)-1)^(-1/theta)/ (1-q))    
    stable <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="stable"),
                      control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    out <- c(out,switch(tail,lwr=(exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)+2*q-1)/q,
                        upr=exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)/(1-q)))
    invgauss <- emfrail(formula = formula, data = data, distribution=emfrail_dist(dist="pvf"), 
                        control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    out <- c(out,switch(tail,
                        lwr=(exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)+2*q-1)/q,
                        upr=exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)/(1-q)))
    ## Non-parametric estimate
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    Names <- c("Estimate", "Gamma", "Positive stable", "Inverse Gaussian")
    id <- rep(1:length(x),2)
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    if(method == "dabrowska"){
        haz <- biHazards(formula, data)
        H <- (haz$lambda10 * haz$lambda01 - haz$lambda11) / ((1 - haz$lambda10) * (1 - haz$lambda01))
        H[is.nan(H)] <- 0
        KMx <- KaplanMeier(x,xstatus)
        KMy <- KaplanMeier(y,ystatus)
        Fx <- 1 - KMx
        Fy <- 1 - KMy
        indx <- match(1, Fx >= q)
        indy <- match(1, Fy >= q)
        qx <- xuni[indx]
        qy <- yuni[indy]
        if(tail == "lwr"){
            prob <- KMx[indx] * KMy[indy] * prod(1 - H[xuni <= qx, yuni <= qy])
            out <- c((2*q - 1 + prob) / q, out)
        }
        else{
            out <- c((KMx[indx] * KMy[indy] * prod(1 - H[xuni <= qx, yuni <= qy])) / (1 - q), out)
        }
    }
    else{
        KMy <- KaplanMeier(y, ystatus)
        FF <- 1 - KMy
        qq <- yuni[match(1, FF >= q)]
        ifelse(tail == "lwr", indy <- (y <= qq & ystatus == 1), indy <- y > qq) 
        xx <- x[indy]
        xxstatus <- xstatus[indy]
        KMx <- KaplanMeier(xx, xxstatus)
        Fx <- 1 - KMx
        ifelse(tail == "lwr", out <- c(Fx[match(1, xuni >= qq) - 1], out), out <- c(KMx[match(1, xuni >= qq)], out))
    }
    if(out[1] < 0 | out[1] > 1){
        stop("q is too close to either 0 or 1")
    }
    ans <- data.frame(Distribution = Names, TailDependence = out)
    colnames(ans) <- c("Distribution", "Tail dependence")
    ans
}

    
#' Bootstrapped confidence interval for "tail dependence"
#'
#' @title Bootstrapped confidence interval for "tail dependence"
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param q Quantile to estimate "tail-dependence" for
#' @param method How to estimate survival probabilities. Available are 'dabrowska' and 'fast'
#' @param tail Tail to estimate "tail dependence" for
#' @param n Number of bootstraps
#' @param level The confidence level required.
#' @param type Type of confidence interval. Available are "normal", which is equal to estimate plus/minus 1.96 times the bootstrap standard error,
#' and "quantile", which is the empirical quantiles of the bootstrapped estimates. 
#' @return CI for "tail dependence"
#' @seealso tailDep
#' @import graphics
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDepCI <- function(formula, data, q, method = "fast", tail="lwr", n = 1000, level = .95, type = "normal"){
    if(!(type %in% c("normal", "quantile")))
        stop("'type' has to be either 'normal' or 'quantile'")
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    data <- data.frame(time = c(x,y), status = c(xstatus,ystatus), cluster = rep(1:length(x),2))
    sims <- numeric(n)
    for(i in 1:n){
        ind <- sample(1:length(x), replace = TRUE)
        d <- data.frame(time = c(x[ind], y[ind]), status = c(xstatus[ind], ystatus[ind]),
                        id = rep(1:length(x),2))
        sims[i] <- tailDep(Surv(time,status) ~ cluster(id), data = d, q, method = method, tail=tail)
    }
    d <- data.frame(time = c(x, y), status = c(xstatus, ystatus), id = rep(1:length(x),2))
    tmp <- tailDep(Surv(time,status) ~ cluster(id), data = d, q, method = method, tail=tail)
    a <- (1-level)/2
    a <- c(a,1-a)
    ## pct <- format.perc(a, 3)
    pct <- paste(a*100, "%")
    ci <- array(NA_real_, dim = c(1L, 2L), dimnames = list("Tail dependence", pct))  
    if(type=="normal"){
        ci[] <- tmp +  qnorm(a) * sd(sims)
    }
    else{
        ci[] <- quantile(a, sims)
    }
    ci
}

#' Function that plots the "tail-dependence" as a function of the chosen quantile
#'
#' @title Function that plots the "tail-dependence" as a function of the chosen quantile.
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data a data.frame containing the variables in the model.
#' @param q Quantiles to estimate "tail-dependence" for.
#' @param tail Tail to estimate "tail dependence" for ("lwr" or "upr").
#' @param method What estimator to use. Can either be "fast" or "dabrowska".
#' @return Plot of "tail-dependence" as a function of the quantile; estimated and implied by frailty models.
#' @seealso tailDepCI tailDependence
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDependencePlot <- function(formula, data, q = seq(.1,.2,.005), tail = "lwr",
                               method = "fast"){
    if(tail == "upr" & all(q == seq(.1,.2,.005))) q <- 1-q
    out <- matrix(NA, nrow = length(q), ncol = 4)
    colnames(out) <- c("Estimate", "Gamma", "Positive stable", "Inverse Gaussian")
    out[,1] <- sapply(q, function(x) tailDep(formula, data, x, tail, method))
    d <- uniTrans(formula, data)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus; id <- rep(1:length(x),2)
    gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
    theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
    out[,2] <- switch(tail,lwr=((2*(1-q)^(-theta) - 1)^(-1/theta) + 2*q - 1) / q,
                      upr=(2*(1-q)^(-theta)-1)^(-1/theta)/ (1-q))
    stable <- emfrail(formula, data, distribution=emfrail_dist(dist="stable"), control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    out[,3] <- switch(tail,lwr=(exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)+2*q-1)/q,
                      upr=exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)/(1-q))
    invgauss <- emfrail(formula, data, distribution=emfrail_dist(dist="pvf"), control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    out[,4] <- switch(tail, lwr=(exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)+2*q-1)/q,
                      upr=exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)/(1-q))
    out <- reshape2::melt(out)
    out <- cbind(out, q)
    names(out[,"Var2"]) <- "Distribution"
    ggplot(mapping = aes(x = .data$q, y = .data$value, group = .data$Var2, colour = .data$Var2), data = out) +
        geom_line() + theme_bw() + xlab("Quantile") + ylab("Tail dependence")
}

tailDep <- function(formula, data, q, tail = "lwr", method = "fast"){
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    Names <- c("Estimate", "Gamma", "Positive stable", "Inverse Gaussian")
    id <- rep(1:length(x),2)
    if(!(method %in% c("dabrowska","fast"))){
        stop("method has to be either 'dabrowska' or 'fast'")
    }
    if(!(tail %in% c("lwr","upr"))){
        stop("tail has to be either 'lwr' or 'upr'")
    }
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    if(method == "dabrowska"){
        haz <- biHazards(formula, data)
        H <- (haz$lambda10 * haz$lambda01 - haz$lambda11) / ((1 - haz$lambda10) * (1 - haz$lambda01))
        H[is.nan(H)] <- 0
        KMx <- KaplanMeier(x,xstatus)
        KMy <- KaplanMeier(y,ystatus)
        Fx <- 1 - KMx
        Fy <- 1 - KMy
        indx <- match(1, Fx >= q)
        indy <- match(1, Fy >= q)
        qx <- xuni[indx]
        qy <- yuni[indy]
        if(tail == "lwr"){
            prob <- KMx[indx] * KMy[indy] * prod(1 - H[xuni <= qx, yuni <= qy])
            out <- (2*q - 1 + prob) / q
        }
        else{
            out <- (KMx[indx] * KMy[indy] * prod(1 - H[xuni <= qx, yuni <= qy])) / (1 - q)
        }
    }
    else{
        KMy <- KaplanMeier(y, ystatus)
        FF <- 1 - KMy
        qq <- yuni[match(1, FF >= q)]
        ifelse(tail == "lwr", indy <- (y <= qq & ystatus == 1), indy <- y > qq) 
        xx <- x[indy]
        xxstatus <- xstatus[indy]
        KMx <- KaplanMeier(xx, xxstatus)
        Fx <- 1 - KMx
        ifelse(tail == "lwr", out <- Fx[match(1, xuni >= qq) - 1], out <- KMx[match(1, xuni >= qq)])
    }
    if(out < 0 | out > 1){
        stop("q is too close to either 0 or 1")
    }
    out
}
