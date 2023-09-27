library(biSurv)
library(survival)
data(diabetic)
data(kidney)
test1 <- independenceTest(Surv(time,status)~cluster(id), data = diabetic)
test1
test2 <- independenceTest(Surv(time,status)~cluster(id), data = kidney)
test2

tmpdatarep <- data.frame(time = rep(rexp(100), 4), id = rep(1:200,2), status = 1)
independenceTest(Surv(time,status)~cluster(id), data = tmpdatarep)
uni <- uniTrans(Surv(time,status)~cluster(id), data = tmpdatarep)
