#' Dabrowska estimator of the bivariate survival function
#'
#' @title Dabrowska estimator of the bivariate survival function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function.
#' The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @return an object of class 'dabrowska' containing the dabrowska estimate of the
#' bivariate survival function
#' @importFrom Rfast sort_unique
#' @examples
#' data(kidney)
#' dabrowska(Surv(time,status) ~ age + sex + disease + cluster(id), data = kidney)
#' @seealso print.dabrowska plot.dabrowska CHR
#' @examples
#' library(survival)
#' data("diabetic")
#' s <- dabrowska(Surv(time,status) ~ cluster(id), data = diabetic)
#' s
#' plot(s)
#' @details the Dabrowska estimator of the bivariate survival function is estimated
#' and returned along with the independence estimator. This function may be interesting
#' in and of itself or because of its use in determining goodness of fit internally in the \code{CHR}
#' function or elsewhere. 
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
dabrowska <- function(formula, data){
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
        t01 <- NA
        t02 <- NA
        medsurv <- est[length(xuni),length(yuni)]
    }
    medCon <- 4 * medsurv - 1
    out <- list(surv = est, timex = xuni, timey = yuni,
                indep = indep, n = nrow(data), n.events = sum(xstatus+ystatus),
                medsurv = medsurv, medCon = medCon, median1 = t01, median2 = t02)
    class(out) <- "biSurv"
    out
}

#' Print a short summary of dabrowska estimate
#' @title Print a short summary of dabrowska estimate
#' @param x an object of class "dabrowska"
#' @param digits minimal number of _significant_ digits, see \code{print.default}.
#' @param ... further arguments for \code{ggplot}.
#' @return number of observations, number of events, medians for the marginals,
#' median concordance, and joint probability of both marginals being greater than respective medians
#' @references Hougaard, Philip. (2000). Analysis of Multivariate Survival Data.
#' @references Dabrowska, Dorota M. "Kaplan-Meier estimate on the plane." The Annals of Statistics (1988): 1475-1489.
#' @seealso dabrowska plot.dabrowska
#' @examples
#' library(survival)
#' data("diabetic")
#' s <- dabrowska(Surv(time,status)~cluster(id), data = diabetic)
#' s
#' plot(s)
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
print.biSurv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    print(data.frame(n = x$n, events = x$n.events, median1 = x$median1, median2 = x$median2,
                     concordance = x$medCon, survival = x$medsurv))
}

#' Contour plot of the Dabrowska estimate of the bivariate survival function
#'
#' @title Contour plot of the Dabrowska estimate of the bivariate survival function
#' @param x an object of class \code{dabrowska}.
#' @param ... further arguments for \code{ggplot}.
#' @return contour plot of estimate of bivariate survival function
#' along with contour plot of independence estimate of bivariate survival function.  
#' @references Hougaard, Philip. (2000). Analysis of Multivariate Survival Data.
#' @references Dabrowska, Dorota M. "Kaplan-Meier estimate on the plane." The Annals of Statistics (1988): 1475-1489.
#' @seealso dabrowska print.dabrowska
#' @examples
#' library(survival)
#' data("diabetic")
#' s <- dabrowska(Surv(time,status)~cluster(id), data = diabetic)
#' plot(s)
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
plot.biSurv <- function(x, ...){
    par(mfrow=c(1,2))
    contour(x = x$timex, y = x$timey, z = x$surv, xlab = "", ylab = "", main = "Estimate")
    contour(x = x$timex, y = x$timey, z = x$indep, xlab = "", ylab = "", main = "Independence")
}

    
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

KaplanMeier <- function(time, status){
    t <- sort_unique(time)
    Y <- sapply(t, function(x) sum(time>=x))
    dN <- sapply(t, function(x) sum(time==x & status == 1))
    cumprod(1-dN/Y)
}

NPMLE <- function(formula, data, maxIt = 100){
    tol <- .Machine$double.eps
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- c(sort_unique(x), max(x)+1)
    yuni <- c(sort_unique(y), max(x)+1)
    n <- length(x)
    n1 <- length(xuni)
    n2 <- length(yuni)
    p1 <- p3 <- matrix(0, nrow = n1, ncol = n2)
    
### Step 1 ###
    for(i in 1:n){
        if(xstatus[i] == 1 & ystatus[i] == 1){
            p3[xuni == x[i], yuni == y[i]] <- p3[xuni == x[i], yuni == y[i]] + 1 / n
        }
        else if(xstatus[i] == 1 & ystatus[i] == 0){
            p1[xuni == x[i], yuni > y[i]] <- p1[xuni == x[i], yuni > y[i]] + 
                1 / (n * sum(yuni > y[i]))
        }
        else if(xstatus[i] == 0 & ystatus[i] == 1){
            p1[xuni > x[i], yuni == y[i]] <- p1[xuni > x[i], yuni == y[i]] + 
                1 / (n * sum(xuni > x[i]))
        }
        else{
            p1[xuni > x[i], yuni > y[i]] <- p1[xuni > x[i], yuni > y[i]] +
                1 / (n * sum(xuni > x[i]) * sum(yuni > y[i]))
        }
    }
    p1 <- p1 + p3
    
### Step 2 ###
    error <- 1
    itt <- 0
    indx01 <- which(xstatus == 0 & ystatus == 1)
    indx10 <- which(xstatus == 1 & ystatus == 0)
    indx00 <- which(xstatus == 0 & ystatus == 0)
    while(error > tol & itt < maxIt){
        itt <- itt + 1
        p2 <- p3
        for(i in indx01){
            tmp <- (xuni > x[i])
            p2[tmp, yuni == y[i]] <- p2[tmp, yuni == y[i]] + 1/n *
                p1[tmp, yuni == y[i]] / sum(p1[tmp, yuni == y[i]])
        }
        for(i in indx10){
            tmp <- (yuni > y[i])
            p2[xuni == x[i], tmp] <- p2[xuni == x[i], tmp] + 1/n *
                p1[xuni == x[i], tmp] / sum(p1[xuni == x[i], tmp])                
        }
        for(i in indx00){
                p2[xuni>x[i], yuni>y[i]] <- p2[xuni>x[i], yuni>y[i]] + 1/n * 
                    p1[xuni>x[i], yuni>y[i]] / sum(p1[xuni>x[i], yuni>y[i]])
        }
        error <- max(abs(p2 - p1))
        p1 <- p2
    }

### Last stuff ###
    cat("Number of iterations: ", itt, "\n")
    KMx <- KaplanMeier(x, xstatus)
    KMy <- KaplanMeier(y, ystatus)   
    indep <- KMx %*% t(KMy)
    est <- t(apply(apply(p1[n1:1,n2:1], 2, cumsum), 
                   1, cumsum))[n1:1,n2:1][-1,-1]
    if(any(KMx<.5) & any(KMy<.5)){
        t01 <- match(1, KMx<.5) - 1
        t02 <- match(1, KMy<.5) - 1
        medsurv <- est[t01,t02] 
    }
    else{
        t01 <- NA
        t02 <- NA
        medsurv <- est[length(xuni)-1,length(yuni)-1]
    }
    medCon <- 4 * medsurv - 1
    out <- list(surv = est, timex = xuni[-length(xuni)], timey = yuni[-length(yuni)],
                indep = indep, n = nrow(data), n.events = sum(xstatus+ystatus),
                medsurv = medsurv, medCon = medCon, median1 = t01, median2 = t02)
    class(out) <- "biSurv"
    out
}

pruitt <- function(formula, data, gamma, maxIt = 100){
    tol <- .Machine$double.eps
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- c(sort_unique(x), max(x)+1)
    yuni <- c(sort_unique(y), max(x)+1)
    n <- length(x)
    n1 <- length(xuni)
    n2 <- length(yuni)
    p1 <- p3 <- matrix(0, nrow = n1, ncol = n2)
    noneighbour <- 0
  
### Step 1 ###
    for(i in 1:n){
        if(xstatus[i] == 1 & ystatus[i] == 1){
            p3[xuni == x[i], yuni == y[i]] <- p3[xuni == x[i], yuni == y[i]] + 1 / n
        }
        else if(xstatus[i] == 1 & ystatus[i] == 0){
            tmp <- (x[i] - gamma <= x & x[i] + gamma >= x & 
                    y > y[i] & xstatus == 1)
            if(any(tmp)){
                ytmp <- y[tmp]
                ystatustmp <- ystatus[tmp]
                tmpKM <- KaplanMeier(ytmp,ystatustmp)
                tmpprob <- -diff(c(1, tmpKM))
                p3[xuni == x[i], yuni %in% ytmp] <- p3[xuni == x[i], yuni %in% ytmp] + tmpprob / n
                p3[xuni == x[i], yuni == n2] <- p3[xuni == x[i], yuni == n2] + 1/n*(1-sum(tmpprob))
            }
            else{
                noneighbour <- noneighbour + 1
                p3[xuni == x[i], n2] <- p3[xuni == x[i], n2] + 1/n
            }
        }
        else if(xstatus[i] == 0 & ystatus[i] == 1){
            tmp <- (y[i] - gamma <= y & y[i] + gamma >= y & 
                      x > x[i] & ystatus == 1)
            if(any(tmp)){
                xtmp <- x[y[i] - gamma <= y & y[i] + gamma >= y & 
                          x >= x[i] & ystatus == 1]
                xstatustmp <- xstatus[y[i] - gamma <= y & y[i] + gamma >= y & 
                                      x >= x[i] & ystatus == 1]
                tmpKM <- KaplanMeier(xtmp,xstatustmp)
                tmpprob <- -diff(c(1, tmpKM))
                p3[xuni %in% xtmp, yuni == y[i]] <- p3[xuni %in% xtmp, yuni == y[i]] + tmpprob / n
                p3[n1, yuni == y[i]] <- p3[n1, yuni == y[i]] + 1/n*(1-sum(tmpprob))
            }
            else{
                noneighbour <- noneighbour + 1
                p3[n1, yuni == y[i]] <- p3[n1, yuni == y[i]] + 1/n
            }
        }
        else{
                p1[xuni > x[i], yuni > y[i]] <- p1[xuni > x[i], yuni > y[i]] +
                    1 / (n * sum(xuni > x[i]) * sum(yuni > y[i]))
        }
    }
    p1 <- p1 + p3
    cat("Number of observations with too small bandwidth: ", noneighbour, "\n")
    
### Step 2 ###
    error <- 1
    itt <- 0
    inx <- which(xstatus == 0 & ystatus == 0)
    while(error > tol & itt < maxIt){
        itt <- itt + 1
        p2 <- p3
        for(i in inx){
                p2[xuni>x[i], yuni>y[i]] <- p2[xuni>x[i], yuni>y[i]] + 1/n * 
                    p1[xuni>x[i], yuni>y[i]] / sum(p1[xuni>x[i], yuni>y[i]])
        }
        error <- max(abs(p2 - p1))
        p1 <- p2
    }

### Last stuff ###    
    cat("Number of iterations: ", itt, "\n")
    KMx <- KaplanMeier(x, xstatus)
    KMy <- KaplanMeier(y, ystatus)   
    indep <- KMx %*% t(KMy)
    est <- t(apply(apply(p1[n1:1,n2:1], 2, cumsum), 
                   1, cumsum))[n1:1,n2:1][-1,-1]
    if(any(KMx<.5) & any(KMy<.5)){
        t01 <- match(1, KMx<.5) - 1
        t02 <- match(1, KMy<.5) - 1
        medsurv <- est[t01,t02] 
    }
    else{
        t01 <- NA
        t02 <- NA
        medsurv <- est[length(xuni)-1,length(yuni)-1]
    }
    medCon <- 4 * medsurv - 1
    out <- list(surv = est, timex = xuni[-length(xuni)], timey = yuni[-length(yuni)],
                indep = indep, n = nrow(data), n.events = sum(xstatus+ystatus),
                medsurv = medsurv, medCon = medCon, median1 = t01, median2 = t02)
}


#' Non-parametric estimation of the bivariate survival function
#'
#' @title Non-parametric estimation of the bivariate survival function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function.
#' The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param gamma bandwidth of pruitt estimator. Can be ignored for other estimators. 
#' @param maxIt maximum number of iterations.
#' @param method which estimator to use. 'dabrowska' (default), 'pruitt' or NPMLE
#' @return matrix with estimate of bivariate survival function
#' @seealso print.biSurv plot.biSurv
#' @details TODO
#' @references TODO
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
biSurv <- function(formula, data, gamma = NULL, maxIt = 100, method = "dabrowska"){
    Call <- match.call()
    if(method == "dabrowska"){
        hm <- dabrowska(formula, data)
        out <- list(call = Call, surv = hm$surv, timex = hm$timex, timey = hm$timey,
                    indep = hm$indep, n = nrow(data), n.events = hm$n.events,
                    medsurv = hm$medsurv, medCon = hm$medCon, median1 = hm$median1,
                    median2 = hm$median2)
        class(out) <- "biSurv"
        out
    }
    else if(method == "NPMLE"){
        hm <- NPMLE(formula, data, maxIt = maxIt)
        out <- list(call = Call, surv = hm$surv, timex = hm$timex, timey = hm$timey,
                    indep = hm$indep, n = nrow(data), n.events = hm$n.events,
                    medsurv = hm$medsurv, medCon = hm$medCon, median1 = hm$median1,
                    median2 = hm$median2)
        class(out) <- "biSurv"
        out
    }
    else if(method == "pruitt"){
        if(is.null(gamma)){
            stop("you have to specify a bandwidth (gamma)")
        }
        else{
            hm <- pruitt(formula, data, gamma, maxIt = maxIt)
            out <- list(call = Call, surv = hm$surv, timex = hm$timex, timey = hm$timey,
                        indep = hm$indep, n = nrow(data), n.events = hm$n.events,
                        medsurv = hm$medsurv, medCon = hm$medCon, median1 = hm$median1,
                        median2 = hm$median2)
            class(out) <- "biSurv"
            out
        }        
    }
    else{
        stop("method has to be either 'dabrowska', 'pruitt' or 'NPMLE'")
    }
}


