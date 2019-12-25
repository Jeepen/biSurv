#' Estimators of "tail dependence" for bivariate survival data
#'
#' @title Estimators of "tail-dependence" for bivariate survival data
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame in which to interpret the variables named in the
#'          \code{formula} argument.
#' @param q Quantile to estimate "tail-dependence" for
#' @param tail Tail to estimate "tail dependence" for
#' @param method What estimator to use
#' @return Estimate of "tail dependence"
#' @seealso tailDepCI
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDependence <- function(formula, data,q,tail = "lwr", method = "dabrowska"){
    Names <- c("Estimate", "Gamma", "Positive stable", "Inverse Gaussian")
    d <- uniTrans(formula, data)
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
        haz <- biHazards(x,y,xstatus,ystatus)
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
    gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
    theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
    out <- c(out,switch(tail,lwr=((2*(1-q)^(-theta) - 1)^(-1/theta) + 2*q - 1) / q,
                        upr=(2*(1-q)^(-theta)-1)^(-1/theta)/ (1-q)))
    
    stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                      control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
    out <- c(out,switch(tail,lwr=(exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)+2*q-1)/q,
                        upr=exp(-(2*(-log(1-q))^(1/alpha1))^alpha1)/(1-q)))
    invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                        control = emfrail_control(se = F, lik_ci = F, ca_test = F))
    alpha2 <- exp(invgauss$logtheta)
    out <- c(out,switch(tail,
                        lwr=(exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)+2*q-1)/q,
                        upr=exp(alpha2-(alpha2^2+2*log(1-q)*(log(1-q)-2*alpha2))^.5)/(1-q)))
    data.frame(Distribution = Names, TailDependence = out)
}

