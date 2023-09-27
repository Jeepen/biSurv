postmean <- function(formula, data = NULL){
    ## Sorting and stuff incl Cox model
    mf <- match.call()
    Terms <- terms(formula, "cluster", data = data)
    mf$formula <- Terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    m <- grep("cluster", names(mf))
    if(length(m) != 0) cluster <- mf[,m] else stop("formula needs a cluster term")
    X <- model.matrix(formula, data)
    hm <- grep("cluster", colnames(X))
    clust <- X[,hm]
    X <- scale(X[order(clust),-hm], scale = FALSE)
    Y <- model.extract(mf, "response")
    if (!is.Surv(Y)) 
        stop("Response must be a survival object")
    Y <- Y[order(clust),]
    model <- coxph(Y ~ X)
    lp <- model$linear.predictors
    cumhaz <- survfit(model)$cumhaz
    time <- Y[,1]
    status <- Y[,2]
    ## Needed quantities
    timeuni <- sort_unique(time)
    nuni <- length(timeuni)
    n <- length(time)
    ny <- NY(time, status, timeuni, nuni, n)
    N <- ny$N
    Y <- ny$Y
    dN <- diff(rbind(0,N))
    dNdot <- diff(c(0, apply(N, 1, sum)))
    S0 <- sapply(1:nuni, function(t) sum(Y[t,] * exp(lp)) / (n/2))
    Lambda <- sapply(1:nuni, function(t) sum(dNdot[1:t] / (n/2 * S0[1:t])))
    dLambda <- diff(c(0, Lambda))
    H <- sapply(1:n, function(i) sum(Y[,i] * exp(lp[i]) * dLambda))
    R <- function(theta) sapply(1:(n/2), function(i) exp(theta * H[2*i-1]) + exp(theta * H[2*i]) - 1)
    d <- sapply(1:(n/2), function(i) status[2*i-1] + status[2*i])
    Ndot <- sapply(1:(n/2), function(i) N[,2*i-1] + N[,2*i])
    dNdot <- diff(rbind(0, Ndot))
    ## ell <- function(theta){
    ##     ## R <- R(theta)
    ##     (sum((d == 2) * log(1 + theta)) + theta * sum(N[nuni,] * H) - 
    ##      sum((1/theta + Ndot[nuni,]) * log(R(theta)))) / (n/2)
    ## }
    S <- function(theta){
        sums <- sapply(1:(n/2), function(i) H[2*i-1] * exp(theta * H[2*i-1]) + 
                                            H[2*i] * exp(theta * H[2*i]))
        sum2 <- sapply(1:(n/2), function(i) H[2*i-1]^2 + H[2*i]^2)
        if(theta != 0){
            R <- R(theta)
            out <- (sum((d == 2) / (1 + theta)) + sum(log(R)) / theta^2 - 
                    sum((1/theta + Ndot[nuni,]) / R * sums) + sum(N[nuni,] * H)) / (n/2)
        }
        else{
            out <- (sum(d == 2) - sum(Ndot[nuni,] * sums) + 
                    0.5*sum((sums^2 - sum2)) + sum(N[nuni,] * H)) / (n/2)
        } 
        out
    }
    ## theta <- newton(S, x0 = .5)$root
    ## sapply(seq(-2,2,length.out=50), function(x) S(x))
    ## if(S(1e-10)<0){theta <- 0} else{theta <- uniroot(S, interval = c(1e-10, 10))$root}
    theta <- uniroot(S, interval = c(0, 1), extendInt = "yes")$root
    II <- - (S(theta+.0001) - S(theta-.0001)) / .0002
    RR <- R(theta)
    eps <- sapply(1:(n/2), function(i) (d[i] == 2) / (1 + theta) + N[nuni,2*i-1] * H[2*i-1] + 
                                       N[nuni,2*i] * H[2*i] - (1/theta + Ndot[nuni,i]) / RR[i] * 
                                       (H[2*i-1] * exp(theta * H[2*i-1]) + H[2*i] * exp(theta * H[2*i])) + 
                                       log(RR[i]) / theta^2)
    pihelp <- pihelper(lp, H, RR, theta, Ndot, Y, N, nuni, n)
    pi <- apply(pihelp, 1, mean)



    Fhelp <- sapply(1:(n/2), function(i) X[2*i-1,] * H[2*i-1] *
            (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i-1]) *
             exp(theta * H[2*i-1]) + exp(theta * H[2*i-1]) / (theta * RR[i]) +
             N[nuni,2*i-1] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i-1]) /
             RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                        H[2*i] * exp(theta * H[2*i]))) +
            X[2*i,] * H[2*i] *
            (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i]) *
             exp(theta * H[2*i]) + exp(theta * H[2*i]) / (theta * RR[i]) +
             N[nuni,2*i] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i]) /
             RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                        H[2*i] * exp(theta * H[2*i]))))
    ifelse(is.null(dim(Fhelp)), FF <- 0, FF <- apply(Fhelp, 1, sum))
    M <- martin(time, status, timeuni, cumhaz, n, nuni)
    dM <- diff(rbind(0, M))

    SEV <- SEVcpp(lp, Y, X, nuni, n, ncol(X))
    SEV <- function(beta, lp){
        S0 <- sapply(1:nuni, function(t) sum(Y[t,] * exp(lp)) / (n/2))
        S1 <- matrix(NA, nrow = nuni, ncol = ncol(X))
        for(t in 1:nuni){
            for(i in 1:ncol(X)){
                tmp <- apply(Y[t,] * exp(lp) *  
            }
        }
        S1 <- sapply(1:nuni, function(t) sum(Y[t,] * exp(lp)) * X[i,] / (n/2), )                     
        S0 <- S1 <- S2 <- numeric(nuni)
        for(t in 1:nuni){
            S0[t] <- sum(Y[t,] * exp(beta * birthYear)) / (n/2)
            S1[t] <- sum(Y[t,] * exp(beta * birthYear) * birthYear) / (n/2)
            S2[t] <- sum(Y[t,] * exp(beta * birthYear) * birthYear^2) / (n/2)
        }
        E <- S1 / S0
        V <- S2 / S0 - E^2
        list(S0 = S0, S1 = S1, S2 = S2, E = E, V = V)
    }

    
    ## EE <- SEV(beta)$E
    ## w <- numeric(n / 2)
    ## for(i in 1:(n/2)){
    ##     w[i] <- sum((birthYear[2*i-1] - EE) * dM[,2*i-1]) +
    ##         sum((birthYear[2*i] - EE) * dM[,2*i])
    ## } 
    ## sigma1 <- sum(SEV(beta)$V * SEV(beta)$S0 * dLambda)
    ## r <- numeric(nuni)
    ## for(t in 1:nuni) r[t] <- -sum(EE[1:t] * dLambda[1:t])
    ## S0 <- SEV(beta)$S0
    ## gaffel <- matrix(NA, nrow = nuni, ncol = n / 2)
    ## for(t in 1:nuni){
    ##     for(i in 1:(n/2)){
    ##         gaffel[t,i] <- sum(dM[1:t,2*i-1] / S0[1:t]) + 
    ##             sum(dM[1:t,2*i] / S0[1:t]) + r[t] / sigma1 * w[i]
    ##     }
    ## }
    ## dgaffel <- diff(rbind(0, gaffel))

    
    
    list(theta, II, sapply(seq(-.1,.1,by=.01), function(x) S(x)), S(theta), SEV)
}


library(Rcpp)
sourceCpp("~/Dropbox/RPackage/WorkInProgress/postmean.cpp")
sourceCpp("~/Dropbox/RPackage/biSurv/src/martingale.cpp")
source("~/Dropbox/RPackage/biSurv/R/simFrail.R")
library(survival)
library(Rfast)
n <- 100
d <- simFrail(n = n, par = 0.0, cov1 = rep(c(1, 0), n/2), cov2 = rep(c(1, 0), n/2), beta = 2)
dat <- with(d, data.frame(time = c(x,y), status = 1, sex = c(cov1,cov2), id = rep(1:n,2),
                          covari = rnorm(2*n)))
postmean(Surv(time,status) ~ cluster(id),data=dat)
postmean(Surv(time,status) ~ sex + cluster(id),data=dat)
postmean(Surv(time,status) ~ sex + covari + cluster(id),data=dat)
help <- postmean(Surv(time,status)~sex + cluster(id),data=dat)
summary(coxph(Surv(time,status)~sex,data=dat))
coxph(Surv(time,status)~sex,data=dat)$linear.predictors
sourceCpp('gg2.cpp')
hmm <- emfrail(Surv(time,status)~ sex + cluster(id),data=dat)
hmm2 <- cox.aalen(Surv(time,status)~ prop(sex) + cluster(id), data=dat, resample.iid=1)
two.stage(hmm2, data = dat)

## S and E functions
SEV <- function(beta){
  S0 <- S1 <- S2 <- numeric(nuni)
  for(t in 1:nuni){
    S0[t] <- sum(Y[t,] * exp(beta * birthYear)) / (n/2)
    S1[t] <- sum(Y[t,] * exp(beta * birthYear) * birthYear) / (n/2)
    S2[t] <- sum(Y[t,] * exp(beta * birthYear) * birthYear^2) / (n/2)
  }
  E <- S1 / S0
  V <- S2 / S0 - E^2
  list(S0 = S0, S1 = S1, S2 = S2, E = E, V = V)
}
U <- function(beta){
  x <- 0
  E <- SEV(beta)$E
  for(i in 1:n){
    x <- x + sum((birthYear[i] - E) * dN[,i])
  }
  x
}
beta <- uniroot(U, interval = c(-1,1))$root

## Lambda, H
dNdot <- diff(c(0, apply(N, 1, sum)))
Lambda <- numeric(nuni)
S0 <- SEV(beta)$S0
for(t in 1:nuni) Lambda[t] <- sum(dNdot[1:t] / (n/2 * S0[1:t]))
dLambda <- diff(c(0, Lambda))
H <- numeric(n)
for(i in 1:n) H[i] <- sum(Y[,i] * exp(beta * birthYear[i]) * dLambda)
R <- function(theta){
  x <- numeric(n / 2)
  for(i in 1:(n/2)){
    x[i] <- exp(theta * H[2*i-1]) + exp(theta * H[2*i]) - 1
  }
  x
}

## loglikelihood
d <- numeric(n / 2)
for(i in 1:(n/2)) d[i] <- status[2*i-1] + status[2*i]
Ndot <- matrix(NA, nrow = nuni, ncol = n / 2)
# for(t in 1:nuni){
#   for(i in 1:(n/2)){
#     Ndot[t,i] <- N[t,2*i-1] + N[t,2*i]
#   }
# }
for(i in 1:(n/2)) Ndot[,i] <- N[,2*i-1] + N[,2*i]
dNdot <- diff(rbind(0, Ndot))
ell <- function(theta){
  R <- R(theta)
  (sum((d == 2) * log(1 + theta)) + theta * sum(N[nuni,] * H) - 
    sum((1/theta + Ndot[nuni,]) * log(R))) / (n/2)
}
S <- function(theta){
  R <- R(theta)
  sums <- numeric(n / 2)
  for(i in 1:(n/2)){
    sums[i] <- H[2*i-1] * exp(theta * H[2*i-1]) + 
      H[2*i] * exp(theta * H[2*i])
  } 
  (sum((d == 2) / (1 + theta)) + sum(log(R)) / theta^2 - 
    sum((1/theta + Ndot[nuni,]) / R * sums) + sum(N[nuni,] * H)) / (n/2)
}
theta <- uniroot(S, interval = c(1e-3,2))$root
II <- - (S(theta+.0001) - S(theta-.0001)) / .0002

## eps
RR <- R(theta)
eps <- numeric(n / 2)
for(i in 1:(n/2)){
  eps[i] <- (d[i] == 2) / (1 + theta) + N[nuni,2*i-1] * H[2*i-1] + 
    N[nuni,2*i] * H[2*i] - (1/theta + Ndot[nuni,i]) / RR[i] * 
    (H[2*i-1] * exp(theta * H[2*i-1]) + H[2*i] * exp(theta * H[2*i])) + 
    log(RR[i]) / theta^2
}

## pi
pihelp <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    pihelp[t,i] <- exp(beta * birthYear[2*i-1]) * Y[t,2*i-1] * 
      (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i-1]) *
         exp(theta * H[2*i-1]) + exp(theta * H[2*i-1]) / (theta * RR[i]) +
         N[nuni,2*i-1] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i-1]) /
         RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                      H[2*i] * exp(theta * H[2*i]))) +
      exp(beta * birthYear[2*i]) * Y[t,2*i] * 
      (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i]) *
         exp(theta * H[2*i]) + exp(theta * H[2*i]) / (theta * RR[i]) +
         N[nuni,2*i] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i]) /
         RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                      H[2*i] * exp(theta * H[2*i])))
  }
}
pi <- apply(pihelp, 1, mean)

## F
Fhelp <- numeric(n / 2)
for(i in 1:(n/2)){
  Fhelp[i] <- birthYear[2*i-1] * H[2*i-1] *
    (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i-1]) *
       exp(theta * H[2*i-1]) + exp(theta * H[2*i-1]) / (theta * RR[i]) +
       N[nuni,2*i-1] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i-1]) /
       RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                    H[2*i] * exp(theta * H[2*i]))) +
    birthYear[2*i] * H[2*i] *
    (-(1/theta + Ndot[nuni,i]) / RR[i] * (1 + theta * H[2*i]) *
       exp(theta * H[2*i]) + exp(theta * H[2*i]) / (theta * RR[i]) +
       N[nuni,2*i] + (1 + theta * Ndot[nuni,i]) * exp(theta * H[2*i]) /
       RR[i]^2 * (H[2*i-1] * exp(theta * H[2*i-1]) + 
                    H[2*i] * exp(theta * H[2*i])))
}
FF <- mean(Fhelp)

## M
M <- matrix(NA, nrow = nuni, ncol = n)
for(t in 1:nuni){
  for(i in 1:n){
    M[t,i] <- N[t,i] - sum(Y[1:t,i] * dLambda[1:t] * exp(beta * birthYear[i]))
  }
}
dM <- diff(rbind(0, M))
EE <- SEV(beta)$E
w <- numeric(n / 2)
for(i in 1:(n/2)){
  w[i] <- sum((birthYear[2*i-1] - EE) * dM[,2*i-1]) +
    sum((birthYear[2*i] - EE) * dM[,2*i])
} 
sigma1 <- sum(SEV(beta)$V * SEV(beta)$S0 * dLambda)
r <- numeric(nuni)
for(t in 1:nuni) r[t] <- -sum(EE[1:t] * dLambda[1:t])
S0 <- SEV(beta)$S0
gaffel <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    gaffel[t,i] <- sum(dM[1:t,2*i-1] / S0[1:t]) + 
      sum(dM[1:t,2*i] / S0[1:t]) + r[t] / sigma1 * w[i]
  }
}
dgaffel <- diff(rbind(0, gaffel))

## dansetrold
dansetrold <- numeric(n / 2)
for(i in 1:(n/2)){
  dansetrold[i] <- eps[i] + sum(pi * dgaffel[,i]) + FF / sigma1 * w[i]
}

## xi and test statistic
Lambdamins <- matrix(NA, nrow = nuni, ncol = n)
for(t in 1:nuni){
  for(i in 1:n){
    Lambdamins[t,i] <- sum(dLambda[1:t] * Y[1:t,i])
  }
}
RRR <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    RRR[t,i] <- exp(theta*Lambdamins[t,2*i-1]*exp(beta * birthYear[2*i-1])) +
                 exp(theta*Lambdamins[t,2*i]*exp(beta * birthYear[2*i])) - 1
  }
}
xi <- (1 + theta * Ndot) / RRR
W <- apply(xi - 1, 1, sum) / sqrt(n / 2)
supW <- max(abs(W))
#qplot(timeuni, W, geom = "step")

## f, h
fhelp <- hhelp <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    fhelp[t,i] <- 1/RRR[t,i] * (Ndot[t,i] - xi[t,i] * 
        (exp(theta * Lambdamins[t,2*i-1] * exp(beta * birthYear[2*i-1])) * 
           Lambdamins[t,2*i-1] * exp(beta * birthYear[2*i-1]) + 
           exp(theta * Lambdamins[t,2*i]*exp(beta*birthYear[2*i])) *
           Lambdamins[t,2*i] * exp(beta * birthYear[2*i])))
    hhelp[t,i] <- theta * xi[t,i] / RRR[t,i] * 
      (exp(theta * Lambdamins[t,2*i-1] * exp(beta * birthYear[2*i-1])) * 
         Lambdamins[t,2*i-1]*exp(beta*birthYear[2*i-1])*birthYear[2*i-1] +
         exp(theta * Lambdamins[t,2*i] * exp(beta * birthYear[2*i])) * 
         Lambdamins[t,2*i]*exp(beta*birthYear[2*i])*birthYear[2*i])
  }
}
f <- apply(fhelp, 1, mean)
h <- -apply(hhelp, 1, mean)

## g
ghelper <- gghelper(Lambdamins, time, birthYear, RRR, xi, Y, theta, beta,
                    nuni, n / 2)
g <- matrix(NA, nuni, nuni)
for(u in 1:nuni){
  for(t in 1:nuni){
    g[u,t] <- -sum(ghelper[u,t,]) / (n / 2)
  }
}

## Phi
HH <- theta / RRR
Lambdaikl <- matrix(NA, nrow = nuni, ncol = n)
for(i in 1:n) Lambdaikl[,i] <- Lambda * exp(beta * birthYear[i])
lambdaikl <- diff(rbind(0, Lambdaikl))
MM <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    MM[t,i] <- Ndot[t,i] - sum(xi[1:t,i] * Y[1:t,2*i-1] *
      exp(theta * Lambdaikl[1:t,2*i-1]) * lambdaikl[1:t,2*i-1]) -
      sum(xi[1:t,i] * Y[1:t,2*i] *
            exp(theta * Lambdaikl[1:t,2*i]) * lambdaikl[1:t,2*i])
  }
}
dMM <- diff(rbind(0, MM))
Phi <- matrix(NA, nrow = nuni, ncol = n / 2)
for(t in 1:nuni){
  for(i in 1:(n/2)){
    Phi[t,i] <- sum(HH[1:t,i] * dMM[1:t,i]) + f[t] / II * dansetrold[i] +
      sum(g[1:t,t] * dgaffel[1:t,i]) + h[t] / sigma1 * w[i]
  }
}

## Simulations
ZZZ <- matrix(0, nrow = nuni, ncol = n / 2)
sims <- matrix(0, nrow = nuni, ncol = 1000)
for(i in 1:50){
  #print(i)
  ZZZ <- matrix(rep(rnorm(n / 2), nuni), nrow = nuni, ncol = n / 2, byrow = T)
  sims[,i] <- apply(Phi * ZZZ, 1, sum) / sqrt(n / 2)
}
plot(timeuni, W, type = 'l', ylim = c(min(c(sims[,1:50]), W), 
                                      max(c(sims[,1:50]), W)),
     xlab = "Years of Age")
for(i in 1:50) lines(timeuni, sims[,i], col = "gray")
lines(timeuni, W)
# legend("bottomleft", c("sup|W| = 0.03347145", "p-value: 0.5691"))
tests <- apply(abs(sims), 2, max)
pvalues[iii] <- mean(supW < tests)
print(pvalues[iii])
}
summary(pvalues)
qqplot(pvalues, qunif(seq(1:12) / 13))
abline(0,1)
ks.test(pvalues, punif)
sumtest <- sum(qnorm(pvalues)) / sqrt(12)
pnorm(sumtest)
(min(pvalues) > (1 - .95^(1 / 12)))
