#' Estimators of "tail-dependence" for bivariate survival data
#'
#' @title Estimators of "tail-dependence" for bivariate survival data
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param q Quantile to estimate "tail-dependence" for
#' @param tail Tail to estimate tail-dependence for
#' @param method What estimator to use
#' @return Estimate of "tail-dependence"
#' @seealso impTailDep
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
tailDep <- function(x, y, xstatus, ystatus, q, tail = "lwr", method = "dabrowska"){
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
            (2*q - 1 + prob) / q
        },upr={
            qx <- min(KMx$time[Fx >= q])
            qy <- min(KMy$time[Fy >= q])
            (KMx$surv[KMx$time == qx] * KMy$surv[KMy$time == qy] * prod(1 - H[xuni <= qx, yuni <= qy])) / (1 - q)
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
            max(Fx[KMx$time < qq])
        },upr={
            max(KMx$surv[KMx$time > qq])
        })  
    })
}

