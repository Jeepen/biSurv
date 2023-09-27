## Practicalities
rm(list = ls())
library(biSurv)
library(reshape2)
library(ggplot2)
library(data.table)
library(survival)
library(frailtyEM)
data("kidney")

## Test
d <- data.frame(time = rexp(300), status = rbinom(300, 1, .5), id = rep(1:150,2))
independenceTest(Surv(time,status) ~ cluster(id), data = d)
tauCens(Surv(time,status) ~ cluster(id), data = d)
(tmp <- CHR(Surv(time,status) ~ cluster(id), data = d))
logLikSort(Surv(time,status) ~ cluster(id), data = d)
plot(tmp)
ss <- dabrowska(Surv(time,status) ~ cluster(id), data = d)
copulaPlot(Surv(time,status) ~ cluster(id), data = d)
ss
plot(ss)
tailDependence(Surv(time,status) ~ cluster(id), data = d, q = .2, method = "fast")
tailDependence(Surv(time,status) ~ cluster(id), data = d, q = .2, method = "dabrowska")
tailDepCI(Surv(time,status) ~ cluster(id), data = d, q = .2, method = "fast")
tailDepCI(Surv(time,status) ~ cluster(id), data = d, q = .2, method = "dabrowska")
tailDependencePlot(Surv(time,status) ~ cluster(id), data = d, q = seq(.15,.25,.005), method = "fast")

## Tau
independenceTest(Surv(time,status) ~ age + sex + disease + cluster(id), data = kidney)
independenceTest(Surv(time,status) ~ cluster(id), data = kidney)
tauCens(Surv(time,status) ~ age + sex + disease + cluster(id), data=kidney)
tauCens(Surv(time,status) ~ cluster(id), data=kidney)
tauCens(Surv(time,status) ~ cluster(id), data=kidney, method = "naive")
tmp <- CHR(Surv(time,status) ~ age + sex + disease + cluster(id), data=kidney)
tmp <- CHR(Surv(time,status) ~ cluster(id), data=kidney)
plot(tmp)
logLikSort(Surv(time,status) ~ cluster(id), data=kidney)
uniTrans(Surv(time,status) ~ cluster(id), data = kidney)
tailDependence(Surv(time,status) ~ cluster(id), data = kidney, q = .1)
tailDependencePlot(Surv(time,status) ~ cluster(id), data = kidney, tail = "upr", q = seq(.8,.9,.01))
tailDepCI(Surv(time,status) ~ cluster(id), data = kidney, q = .4, n = 100)

data(diabetes)
independenceTest(Surv(time,status) ~ cluster(id), data = diabetes)
tauCens(Surv(time,status) ~ cluster(id), data = diabetes)
(tmp <- CHR(Surv(time,status) ~ cluster(id), data = diabetes))
plot(tmp)
tailDependencePlot(Surv(time,status) ~ cluster(id), data = diabetes)

## Simulate
d <- simFrail(par = 0, cov1 = rep(c(1,0),50), cov2 = rep(c(1,0),50), beta = 2)
dat <- with(d, data.frame(time = c(x,y), status = 1, sex = c(cov1,cov2), id = rep(1:100,2)))
m0 <- emfrail(Surv(time,status) ~ cluster(id), data = dat)
m1 <- emfrail(Surv(time,status) ~ sex + cluster(id), data = dat)
independenceTest(Surv(time,status) ~ cluster(id), data = dat)
independenceTest(Surv(time,status) ~ sex + cluster(id), data = dat)

tailDependence(Surv(time,status) ~ cluster(id), data = dat, q = .1, method = "fast")
tailDependence(Surv(time,status) ~ cluster(id), data = dat, q = .2, method = "dabrowska")

tauCens(Surv(time,status) ~ cluster(id), data = dat)
tauCens(Surv(time,status) ~ sex + cluster(id), data = dat)

data(diabetes)
m5 <- coxph(Surv(time,status) ~ sex, data = dat)
mm <- cox.aalen(Surv(time,status) ~ prop(sex) + cluster(id), data = dat)
fit <- two.stage(mm, data = dat)
mm <- cox.aalen(Surv(time,status) ~ cluster(id), data = dat)

marg <- cox.aalen(Surv(time,status)~prop(treat)+prop(adult)+
                      cluster(id),data=diabetes,resample.iid=1)
ddd <- postmean(Surv(time,status) ~ cluster(id), data = dat)


sourceCpp("~/Dropbox/RPackage/WorkInProgress/postmean.cpp")
tmp <- with(dat, NY(time, status, sort(unique(time)), length(unique(time)), length(time)))



