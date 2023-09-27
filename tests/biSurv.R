library(biSurv)
library(survival)
library(Rcpp)
data(diabetic)
data(kidney)
library(ggplot2)
library(microbenchmark)
library(data.table)
source("~/Dropbox/RPackage/tests/simfunction.R")

S <- biSurv(Surv(time,status) ~ cluster(id), data = diabetic)
S2 <- biSurv(Surv(time,status) ~ cluster(id), data = diabetic, method = "NPMLE", maxIt = 1e4)
S3 <- biSurv(Surv(time,status) ~ cluster(id), data = diabetic, method = "pruitt", gamma = 10, maxIt = 1e4)
S
S2
S3
plot(S)
plot(S2)
plot(S3)
S$surv[1:5,1:7]
S2$surv[1:5,1:7]
S3$surv[1:5,1:7]
S3$surv[178:184,178:182]
tauCens(Surv(time,status) ~ cluster(id), data = diabetic)

tmp <- uniTrans(Surv(time,status) ~ cluster(id), data = diabetic)

S <- biSurv(Surv(time,status)~cluster(id), data = kidney)
S2 <- biSurv(Surv(time,status)~cluster(id), data = kidney, method = "NPMLE", maxIt = 1000)
S3 <- biSurv(Surv(time,status)~cluster(id), data = kidney, method = "pruitt", gamma = 5, maxIt = 1000)
S
S2
S3
plot(S)
plot(S2)
plot(S3)

tmp <- simFrail(n = 100)
setDT(tmp)
tmp[, C := rexp(200,rate=.5)]
tmp[, status := as.numeric(C>=time)]
tmp[, time := pmin(time, C)]
tmp2 <- uniTrans(Surv(time,status)~cluster(id),data=tmp)
tmp2
surver <- biSurv(Surv(time,status) ~ cluster(id), data = tmp, method = "pruitt", gamma = .1, maxIt = 1000)
surver2 <- biSurv(Surv(time,status) ~ cluster(id), data = tmp, method = "dabrowska")
surver3 <- biSurv(Surv(time,status) ~ cluster(id), data = tmp, method = "NPMLE", maxIt = 10000)
surver2
surver
surver3
plot(surver2)
plot(surver)
plot(surver3)
tauCens(Surv(time,status) ~ cluster(id), data = kidney)
medianConcordance(Surv(time,status) ~ cluster(id), data = kidney)
independenceTest(Surv(time,status) ~ cluster(id), data = kidney, weight = 1)

S1 <- biSurv(Surv(time,status) ~ cluster(id), data = tmp)
plot(S1)
S2 <- biSurv(Surv(time,status) ~ cluster(id), data = tmp, method = "pruitt", gamma = .2)
plot(S2)
S3 <- biSurv(Surv(time,status) ~ cluster(id), data = tmp, method = "NPMLE")
plot(S3)


tmp <- microbenchmark(biSurv(Surv(time,status)~cluster(id), data = diabetic),
                      biSurv(Surv(time,status)~cluster(id), data = diabetic, method = "NPMLE"),
                      biSurv(Surv(time,status)~cluster(id), data = diabetic, method = "pruitt", gamma = 5))




