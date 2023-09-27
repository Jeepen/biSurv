#' EM estimator of bivariate survival function
#'
#' @title EM estimator of bivariate survival function
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function.
#' The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param tol how large is the difference from iteration to iteration allowed to be
#' @param maxIt how many iterations are acceptable
#' @return matrix with estimate of bivariate survival function
#' @seealso dabrowska
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
EMsurv <- function(formula, data, tol = 1e-3, maxIt = 10){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xu <- sort_unique(x)
    yu <- sort_unique(y)
    n <- length(x)
    n1 <- length(xu)
    n2 <- length(yu)
    p1 <- p3 <- matrix(0, nrow = n1, ncol = n2)
    
### Step 1 ###
    for(i in 1:n){
        if(xstatus[i] == 1 & ystatus[i] == 1){
            p3[xu == x[i], yu == y[i]] <- p3[xu == x[i], yu == y[i]] + 1 / n
        }
        else if(xstatus[i] == 1 & ystatus[i] == 0 & any(yu > y[i])){
            p1[xu == x[i], yu > y[i]] <- p1[xu == x[i], yu > y[i]] + 
                1 / (n * sum(yu > y[i]))
        }
        else if(xstatus[i] == 0 & ystatus[i] == 1 & any(xu > x[i])){
            p1[xu > x[i], yu == y[i]] <- p1[xu > x[i], yu == y[i]] + 
                1 / (n * sum(xu > x[i]))
        }
        else if(any(xu>x[i]) & any(yu>y[i])){
            p1[xu > x[i], yu > y[i]] <- p1[xu > x[i], yu > y[i]] +
                1 / (n * sum(xu > x[i]) * sum(yu > y[i]))
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
        print(paste("Iteration: ", itt))
        p2 <- p3
        for(i in indx01){
            tmp <- (xu > x[i])
            if(any(tmp)){
                p2[tmp, yu == y[i]] <- p2[tmp, yu == y[i]] + 1/n *
                    p1[tmp, yu == y[i]] / sum(p1[tmp, yu == y[i]])
            }
        }
        for(i in indx10){
            tmp <- (yu > y[i])
            if(any(tmp)){
                p2[xu == x[i], tmp] <- p2[xu == x[i], tmp] + 1/n *
                    p1[xu == x[i], tmp] / sum(p1[xu == x[i], tmp])                
            }
        }
        for(i in indx00){
            if(any(xu > x[i]) & any(yu > y[i])){
                p2[xu > x[i], yu > y[i]] <- p2[xu > x[i], yu > y[i]] + 1/n *
                    p1[xu > x[i], yu > y[i]] / sum(p1[xu > x[i], yu > y[i]])                
            }
        }
        error <- max(abs(p2 - p1))
        p1 <- p2
        print(error)
    }
    list(p = p1, S = t(apply(apply(p1[n1:1,n2:1], 2, cumsum), 
                             1, cumsum))[n1:1,n2:1])
}

#' Pruitt estimator of bivariate survival function
#'
#' @title Pruitt estimator of bivariate survival function
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param gamma Bandwidth
#' @param tol How large is the difference from iteration to iteration allowed to be
#' @param maxIt How many iterations are acceptable
#' @return Matrix with estimate of bivariate survival function
#' @seealso dabrowska EMsurv
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
pruitt <- function(formula, data, gamma, tol = 1e-3, maxIt = 10){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4)
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xu <- sort_unique(x)
    yu <- sort_unique(y)
    n <- length(x)
    n1 <- length(xu)
    n2 <- length(yu)
    p1 <- p3 <- matrix(0, nrow = n1, ncol = n2)
    noneighbour <- 0
  
### Step 1 ###
    for(i in 1:n){
        if(xstatus[i] == 1 & ystatus[i] == 1){
            p3[xu == x[i], yu == y[i]] <- p3[xu == x[i], yu == y[i]] + 1 / n
        }
        else if(xstatus[i] == 1 & ystatus[i] == 0){
            tmp <- (x[i] - gamma <= x & x[i] + gamma >= x & 
                    y > y[i] & xstatus == 1)
            if(any(tmp)){
                ytmp <- y[tmp]
                ystatustmp <- ystatus[tmp]
                tmpKM <- KaplanMeier(ytmp,ystatustmp)
                tmpprob <- -diff(c(1, tmpKM[-length(tmpKM)],0))
                p3[xu == x[i], yu %in% ytmp] <- p3[xu == x[i], yu %in% ytmp] + tmpprob / n
            }
            else{
                noneighbour <- noneighbour + 1
                warning("no observations in neighbourhood: consider more smoothing (higher gamma)")
                print(noneighbour)
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
                tmpprob <- -diff(c(1, tmpKM[-length(tmpKM)],0))
                p3[xu %in% xtmp, yu == y[i]] <- p3[xu %in% xtmp, yu == y[i]] + tmpprob / n
            }
            else{
                noneighbour <- noneighbour + 1
                warning("no observations in neighbourhood: consider more smoothing (higher gamma)")
                print(noneighbour)
            }
            ## for(w in 1:length(tmpprob)){
            ##     p3[xu == tmpKM$time[w], yu == y[i]] <- 
            ##         p3[xu == tmpKM$time[w], yu == y[i]] + 1 / n * tmpprob[w]
            ## }
        }
        else{
            if(any(xu > x[i]) & any(yu > y[i])){
                p1[xu > x[i], yu > y[i]] <- p1[xu > x[i], yu > y[i]] +
                    1 / (n * sum(xu > x[i]) * sum(yu > y[i]))
            }
        }
    }
    p1 <- p1 + p3
    
### Step 2 ###
    error <- 1
    itt <- 0
    inx <- which(xstatus == 0 & ystatus == 0)
    while(error > tol & itt < maxIt){
        itt <- itt + 1
        print(paste("Iteration: ", itt))
        p2 <- p3
        for(i in inx){
            if(any(xu > x[i]) & any(yu>y[i])){
                p2[xu>x[i], yu>y[i]] <- p2[xu>x[i], yu>y[i]] + 1/n * 
                    p1[xu>x[i], yu>y[i]] / sum(p1[xu>x[i], yu>y[i]])
            }
        }
        error <- max(abs(p2 - p1))
        p1 <- p2
        print(error)
    }
    list(p = p1, 
         S = t(apply(apply(p1[n1:1,n2:1], 2, cumsum), 1, cumsum))[n1:1,n2:1])
}

KaplanMeier <- function(time, status){
    t <- sort_unique(time)
    Y <- sapply(t, function(x) sum(time>=x))
    dN <- sapply(t, function(x) sum(time==x & status == 1))
    cumprod(1-dN/Y)
}
