#' Dabrowska estimator of the bivariate survival function
#'
#' @title Dabrowska estimator of the bivariate survival function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @return Matrix with Dabrowska estimate of the bivariate survival function
#' @importFrom Rfast sort_unique
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
dabrowska <- function(formula, data){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    haz <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
    KMx <- KaplanMeier(x, xstatus)
    KMy <- KaplanMeier(y, ystatus)   
    H <- (haz$lambda10*haz$lambda01-haz$lambda11) / 
        ((1-haz$lambda10)*(1-haz$lambda01))
    H[is.nan(H)] <- 0
    indep <- KMx %*% t(KMy)
    est <- t(apply(apply(1 - H, 2, cumprod), 1, cumprod)) * indep
    if(any(KMx<.5) & any(KMy<.5)){
        t01 <- match(1, KMx<.5) - 1
        t02 <- match(1, KMy<.5) - 1
        medsurv <- est[t01,t02] 
    }
    else{
        t01 <- Inf
        t02 <- Inf
        medsurv <- est[length(xuni),length(yuni)]
    }
    medCon <- 4*medsurv - 1
    out <- list(call = Call, surv = est, timex = xuni, timey = yuni,
                indep = indep, n = nrow(data), n.events = sum(xstatus+ystatus),
                medsurv = medsurv, medCon = medCon, median1 = t01, median2 = t02)
    class(out) <- "dabrowska"
    out
}

#' @export
print.dabrowska <- function(x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor, 
                            signif.stars = getOption("show.signif.stars"), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(data.frame(n = x$n, events = x$n.events, median1 = x$median1, median2 = x$median2, concordance = x$medCon, probability = x$medsurv))
}

#' Contour plot of the Dabrowska estimate of the bivariate survival function
#'
#' @title Contour plot of the Dabrowska estimate of the bivariate survival function
#' @param x An object of class \code{dabrowska}.
#' @param ... Further arguments for \code{ggplot}.
#' @return Contour plot of estimate of bivariate survival function along with contour plot of independence estimate of bivariate survival function.  
#' @seealso dabrowska
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
plot.dabrowska <- function(x, ...){
    par(mfrow=c(1,2))
    contour(x = x$timex, y = x$timey, z = x$surv, xlab = "", ylab = "", main = "Estimate")
    contour(x = x$timex, y = x$timey, z = x$indep, xlab = "", ylab = "", main = "Independence")
}

    
#' @export
biHazards <- function(formula, data = NULL){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    h <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
    h$lambda11[is.nan(h$lambda11)] <- 0
    h$lambda10[is.nan(h$lambda10)] <- 0
    h$lambda01[is.nan(h$lambda01)] <- 0
    h
}


#' @export
KaplanMeier <- function(time, status){
    t <- sort_unique(time)
    Y <- sapply(t, function(x) sum(time>=x))
    dN <- sapply(t, function(x) sum(time==x & status == 1))
    cumprod(1-dN/Y)
}
