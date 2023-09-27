independenceTest <- function(x,y,xstatus,ystatus, weight = "independence"){
    eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
    nu1 <- nrow(eyy)
    nu2 <- ncol(eyy)
    m1 <- survfit(Surv(x,xstatus)~1)
    m2 <- survfit(Surv(y,ystatus)~1)
    if(class(weight) == "matrix"){}
    else if(weight == "independence"){
        weight <- m1$surv %o% m2$surv
    }
    else if(weight == "dabrowska"){
        weight <- dabrowska(x,y,xstatus,ystatus)
    }
    else if(weight == "atrisk"){
        weight <- eyy
    }
    else if(length(weight)==1){
        weight <- matrix(weight, nrow = nu1, ncol = nu2)
    }
    cumhaz1 <- m1$cumhaz
    cumhaz2 <- m2$cumhaz
    n <- length(x)
    M1 <- matrix(NA, nrow = nu1, ncol = n)
    M2 <- matrix(NA, nrow = nu2, ncol = n)
    for(t in 1:nu1){
        for(i in 1:n){
            M1[t,i] <- xstatus[i] * (x[i] <= m1$time[t]) - cumhaz1[t]*(x[i] > m1$time[t]) -
                cumhaz1[which(m1$time == x[i])]*(x[i] <= m1$time[t])
        }
    }
    for(t in 1:nu2){
        for(i in 1:n){
            M2[t,i] <- ystatus[i] * (y[i] <= m2$time[t]) - cumhaz2[t]*(y[i] > m2$time[t]) -
                cumhaz2[which(m2$time == y[i])]*(y[i] <= m2$time[t])
        }
    }
    dM1 <- rbind(M1[1,],diff(M1))
    dM2 <- rbind(M2[1,],diff(M2))
    testStatistic <- 0
    for(i in 1:n){
       testStatistic <- testStatistic + sum((dM1[,i] %o% dM2[,i]) * weight)
    }
    ## Variance estimate
    help <- c(m1$cumhaz[1],diff(m1$cumhaz)) %o% c(m2$cumhaz[1],diff(m2$cumhaz))
    varianceEst <- sum(weight^2 * eyy * help) * n
    ## Present wrapping
    test <- testStatistic / sqrt(varianceEst)
    d <- data.frame(test, sqrt(varianceEst), 2*pnorm(-abs(test)))
    colnames(d) <- c("Test statistic", "SE", "p-value")
    d
}
