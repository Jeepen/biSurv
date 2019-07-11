#' Integrated squared distance (ISD) between empirical and theoretical CHRs
#'
#' @title Integrated squared distance (ISD) between empirical and theoretical CHRs
#' @param x,y Vectors of failure times
#' @param xstatus,ystatus Status indicators for failure times
#' @param n Number of steps for empirical CHR
#' @param without Distributions to ignore
#' @return Data.frame with ISD for different frailty distributions
#' @seealso chrCpp
#' @references Chen, Min-Chi & Bandeen-Roche, Karen. (2005). A Diagnostic for Association in Bivariate Survival Models. Lifetime data analysis. 11. 245-64. 10.1007/s10985-004-0386-8.
#' @useDynLib biSurv
#' @import frailtyEM
#' @import survival
#' @import ggplot2
#' @import reshape2
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
CHR <- function(x,y,xstatus,ystatus,n=5,without=NULL){
    gammaDiff <- stableDiff <- invgaussDiff <- Names <- NULL
    if(length(x)!=length(y)){
        stop("Length of x and y differ")
    }
    if(length(xstatus)!=length(ystatus)){
        stop("Length of xstatus and ystatus differ")
    }
    n0 <- length(x)
    S <- dabrowska(x,y,xstatus,ystatus)
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
    data <- data.frame(S = seq(.01,.97,.01), Empirical = emp)
    time <- c(x,y)
    status <- c(xstatus,ystatus)
    id <- rep(1:length(x),2)
    if(!("gamma" %in% without)){
        gamma <- coxph(Surv(c(x,y),c(xstatus,ystatus)) ~ frailty(id))
        theta <- gamma$history$`frailty(id)`$history[gamma$iter[1],1]
        gammaCHR <- rep(1+theta, 97)
        gammaDiff <- mean((emp-gammaCHR)^2)
        Names <- c(Names, "Gamma")
        data <- cbind(data, Gamma = gammaCHR)
    }
    if(!("posstab" %in% without)){
        stable <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="stable"), data=data.frame(),
                          control = emfrail_control(se = F, lik_ci = F, ca_test = F))
        alpha1 <- exp(stable$logtheta)/(1+exp(stable$logtheta))
        stableCHR <- 1-(1-alpha1)/(alpha1*log(seq(.01,.97,.01)))
        stableDiff <- mean((emp-stableCHR)^2)
        Names <- c(Names, "Positive stable")
        data <- cbind(data, PositiveStable = stableCHR)
    }
    if(!("invgauss" %in% without)){
        invgauss <- emfrail(Surv(c(x,y),c(xstatus,ystatus)) ~ cluster(id), distribution=emfrail_dist(dist="pvf"), data=data.frame(),
                            control = emfrail_control(se = F, lik_ci = F, ca_test = F))
        alpha2 <- exp(invgauss$logtheta)
        invgaussCHR <- 1+1/(alpha2-log(seq(.01,.97,.01)))
        invgaussDiff <- mean((emp-invgaussCHR)^2)
        Names <- c(Names, "Inverse Gaussian")
        data <- cbind(data,InverseGaussian=invgaussCHR)
    }
    print(data.frame(Distribution = Names, ISD = c(gammaDiff, stableDiff, invgaussDiff)))
    melted <- melt(data,id="S")
    colnames(melted)[2] <- "Model"
    ggplot(melted, aes(x=S,y=value,color=Model)) + geom_line() +
        geom_hline(aes(yintercept=1),linetype=2) + xlab("Survival function value") +
        ylab("CHR")
}

