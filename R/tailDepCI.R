#' Bootstrapped confidence interval for "tail dependence"
#'
#' @title Bootstrapped confidence interval for "tail dependence"
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param cluster Cluster variable
#' @param q Quantile to estimate "tail-dependence" for
#' @param method How to estimate survival probabilities. Available are 'dabrowska' and 'fast'
#' @param tail Tail to estimate "tail dependence" for
#' @param n Number of bootstraps
#' @return CI for "tail dependence"
#' @seealso tailDep
#' @import boot
#' @import graphics
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDepCI <- function(formula, data, cluster, q, method = "dabrowska", tail="lwr", n = 1000){
    cluster <- eval(substitute(cluster),data)
    d <- uniTrans(formula, data, cluster)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    data <- data.frame(time = c(x,y), status = c(xstatus,ystatus), cluster = rep(1:length(x),2))
    sims <- numeric(n)
    for(i in 1:n){
        ind <- sample(1:length(x), replace = TRUE)
        d <- data.frame(time = c(x[ind], y[ind]), status = c(xstatus[ind], ystatus[ind]),
                        id = rep(1:length(x),2))
        sims[i] <- tailDep(Surv(time,status) ~ 1, data = d,
                           cluster = id, q, method = method, tail=tail)
    }
    quantile(sims, c(.025,.975))
}

#' @export
tailDep <- function(formula, data, cluster, q, tail = "lwr", method = "dabrowska"){
    Call <- match.call()
    cluster <- eval(substitute(cluster),data)
    Names <- c("Estimate", "Gamma", "Positive stable", "Inverse Gaussian")
    d <- uniTrans(formula, data, cluster)
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    id <- rep(1:length(x),2)
    if(!(method %in% c("dabrowska","fast"))){
        stop("tail has to be either 'lwr' or 'upr'")
    }
    if(!(tail %in% c("lwr","upr"))){
        stop("tail has to be either 'lwr' or 'upr'")
    }
    switch(method, dabrowska={
        xuni <- sort_unique(x)
        yuni <- sort_unique(y)
        haz <- biHazards(formula, data, cluster)
        H <- (haz$lambda10 * haz$lambda01 - haz$lambda11) / ((1 - haz$lambda10) * (1 - haz$lambda01))
        H[is.nan(H)] <- 0
        KMx <- prodlim(Hist(x,xstatus) ~ 1)
        KMy <- prodlim(Hist(y,ystatus) ~ 1)
        Fx <- 1 - KMx$surv
        Fy <- 1 - KMy$surv
        switch(tail, lwr={
            qx <- min(KMx$time[Fx >= q])
            qy <- min(KMy$time[Fy >= q])
            prob <- KMx$surv[KMx$time == qx] * KMy$surv[KMy$time == qy] * prod(1 - H[xuni <= qx, yuni <= qy])
            out <- (2*q - 1 + prob) / q
        },upr={
            qx <- min(KMx$time[Fx >= q])
            qy <- min(KMy$time[Fy >= q])
            out <- (KMx$surv[KMx$time == qx] * KMy$surv[KMy$time == qy] * prod(1 - H[xuni <= qx, yuni <= qy])) / (1 - q)
        })
    },fast={
        KMy <- prodlim(Hist(y, ystatus) ~ 1)
        FF <- 1-KMy$surv
        qq <- min(KMy$time[FF >= q])
        switch(tail, lwr={
            xx <- x[y < qq & ystatus == 1]
            xxstatus <- xstatus[y < qq & ystatus == 1]
        },upr={
            xx <- x[y > qq]
            xxstatus <- xstatus[y > qq]
        })
        KMx <- prodlim(Hist(xx, xxstatus) ~ 1)
        switch(tail, lwr={
            Fx <- 1-KMx$surv
            out <- max(Fx[KMx$time < qq])
        },upr={
            out <- max(KMx$surv[KMx$time > qq])
        })  
    })
    if(out < 0 | out > 1){
        stop("q is too close to either 0 or 1")
    }
    out
}
