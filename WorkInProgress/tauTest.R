setwd("~/Dropbox/RPackage/WorkInProgress/")
library(Rcpp)
sourceCpp("taucpp2.cpp")
source("tauCens2.R")
library(biSurv)
library(microbenchmark)
library(frailtyEM)
library(Rfast)

d <- data.frame(time = rexp(300), status = rbinom(300, 1, .5), id = rep(1:150,2))
tauCens(Surv(time,status) ~ cluster(id), data = d)
tauCens2(Surv(time,status) ~ cluster(id), data = d)
tmp <- microbenchmark(
    tauCens(Surv(time,status) ~ cluster(id), data = d),
    tauCens2(Surv(time,status) ~ cluster(id), data = d)
)
summary(tmp)


a <- matrix(1:9, nrow = 3, ncol = 3)
b <- matrix(1:9, nrow = 3, ncol = 3)
cppFunction("double fisk(arma::mat a, arma::mat b){
arma::mat c = a%b;
double d = sum(pow(a, 2));
return d;
}", depends = "RcppArmadillo")
fisk(a,b)

cppFunction("arma::vec summe(arma::mat x){
arma::vec d = pow(x.as_col(), 2);
return d;
}", depends = "RcppArmadillo")
summe(a)





setwd("~/Dropbox/RPackage/biSurv/")
library(devtools)
check()




















