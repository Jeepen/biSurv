library(biSurv)
library(survival)
library(Rcpp)
setwd("~/Dropbox/RPackage/WorkInProgress/biSurvEstimator/")
sourceCpp("emcpp.cpp")
source("EM.R")
source("pruitt.R")
library(Rfast)
data(diabetic)
data(kidney)
S <- dabrowska(Surv(time,status) ~ cluster(id), data = diabetic)
S2 <- EMsurv(Surv(time,status) ~ cluster(id), data = diabetic)
S3 <- EMsurv2(Surv(time,status) ~ cluster(id), data = diabetic)
S4 <- pruitt(Surv(time,status) ~ cluster(id), data = diabetic, gamma = 15, maxIt = 10)
S$surv[1:5,1:7]
S2$S[1:5,1:7]
S3[1:5,1:7]
S4$S[1:5,1:7]

start1 <- Sys.time()
S2 <- EMsurv(Surv(time,status) ~ cluster(id), data = diabetic, maxIt = 3)
stop1 <- Sys.time()
start2 <- Sys.time()
S3 <- EMsurv2(Surv(time,status) ~ cluster(id), data = diabetic, maxIt = 3)
stop2 <- Sys.time()
difftime(stop1,start1, units = "secs")
difftime(stop2,start2, units = "secs")
as.numeric(difftime(stop1,start1, units = "secs")) /
    as.numeric(difftime(stop2,start2, units = "secs"))


tmp <- uniTrans(Surv(time,status)~cluster(id), data = diabetic)
attach(tmp)
xuni <- sort_unique(x)
yuni <- sort_unique(y)
tmp2 <- EMcpp(x,y,xstatus,ystatus,xuni,yuni)


chr1 <- CHR(Surv(time,status)~cluster(id), data = diabetic)
chr1
autoplot(chr1)
chr1 <- CHR(Surv(time,status)~cluster(id), data = diabetic)
chr1
autoplot(chr1)
l1 <- logLikSort(Surv(time,status)~cluster(id), data = diabetic)
l1

cppFunction("
int fisk(NumericVector x){
return sum(x > x(0) & x>x(1));
}
")

fisk(rnorm(40))


cppFunction("
arma::mat fisk(arma::mat x, arma::mat y){
return x+y;
}
", depends = "RcppArmadillo")
fisk(matrix(1:9,nrow=3), matrix(2:10,nrow=3))



