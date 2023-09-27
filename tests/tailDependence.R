## Practical
setwd("~/Dropbox/RPackage/tests")
library(biSurv)
library(survival)
library(Rcpp)
library(data.table)
library(tidyverse)
source("simfunction.R")
data(diabetic)
data(kidney)

#Gamma
n <- 1000
ests <- numeric(10)
for(i in 1:10){
  cat(i, "\n")

  gamma <- 2
  X1 <- runif(n)
  X2 <- runif(n)
  V <- rgamma(n, 1/gamma, 1/gamma)
  U1 <- (1-log(X1)/V)^(-1/gamma)
  U2 <- (1-log(X2)/V)^(-1/gamma)
  C1 <- runif(n, max = 3)
  C2 <- runif(n, max = 3)
  
  tmp <- data.table(id = rep(1:n,2), time = c(pmin(U1,C1),pmin(U2,C2)), status = c(C1>U1, C2>U2))
  
  ## tmp <- simFrail(n = 500)
  ## setDT(tmp)
  ## tmp[, C := runif(1e3, max = 3)]
  ## tmp[, status := as.numeric(C>=time)]
  ## tmp[, time := pmin(time, C)]
  # form <- Surv(time,status) ~ cluster(id)
  # tailDependence(form, data=tmp, q = .2)[1,2]
  ests[i] <- biSurv:::tailDep(Surv(time,status) ~ cluster(id), data=tmp, q = .1, tail = "lwr", method = "dabrowska")
  cat(mean(ests[1:i]), "\n")
}
CF <- function(x) (2/(1-x)-1)^(-1) +1 - 2*(1-x)
C <- function(x) (2/x-1)^(-1/gamma)
CF(.1) / .1
C(.1) / .1


tailDependence(Surv(time, status) ~ cluster(id), data = kidney, q = .1)
tailDepCI(Surv(time, status) ~ cluster(id), data = kidney, q = .3)


gamma <- 2
X1 <- runif(1000)
X2 <- runif(1000)
V <- rgamma(1000, 1/gamma, 1/gamma)
U1 <- (1-log(X1)/V)^(-1/gamma)
U2 <- (1-log(X2)/V)^(-1/gamma)
cor(U1,U2,method="kendall")

