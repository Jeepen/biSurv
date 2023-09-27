## Practical
setwd("~/Dropbox/RPackage/tests")
library(biSurv)
library(survival)
library(Rcpp)
library(data.table)
library(ggplot2)
source("simfunction.R")
data(diabetic)
data(kidney)


chr1 <- CHR(Surv(time,status)~cluster(id), data = diabetic)
chr1
autoplot(chr1)
l1 <- logLikSort(Surv(time,status)~cluster(id), data = diabetic)
l1

#Gamma
tmp <- simFrail(n = 500)
setDT(tmp)
tmp[, C := rexp(1e3, rate = .5)]
tmp[, status := as.numeric(C>=time)]
tmp[, time := pmin(time, C)]
test <- CHR(Surv(time,status)~cluster(id),data=tmp)
test
autoplot(test)
logLikSort(Surv(time,status) ~ cluster(id), data=tmp)

#Positive stable
tmp <- simFrail(n = 500, dist = "posstab", par = .7)
setDT(tmp)
tmp[, C := rexp(1e3, rate = .5)]
tmp[, status := as.numeric(C>=time)]
tmp[, time := pmin(time, C)]
test <- CHR(Surv(time,status) ~ cluster(id), data = tmp)
test
autoplot(test)
logLikSort(Surv(time,status) ~ cluster(id), data=tmp)



