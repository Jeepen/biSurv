library(survival)
library(biSurv)
library(Rcpp)
sourceCpp("~/Dropbox/RPackage/biSurv/src/eyy.cpp")

nsim <- 1e2
pvalues <- numeric(nsim)
n <- 1e2

for(i in 1:nsim){
    print(i)
    x <- rexp(n)
    y <- rexp(n)
    cxs <- rexp(n,rate=.5)
    cys <- rexp(n,rate=.5)
    xs <- (x<=cxs)
    ys <- (y<=cys)
    x <- pmin(x,cxs)
    y <- pmin(y,cys)
    ## modelx <- survfit(Surv(x,xs)~1,ctype=1)
    ## resx <- xs - modelx$cumhaz[order(x)]
    ## modely <- survfit(Surv(y,ys)~1,ctype=1)
    ## resy <- ys - modely$cumhaz[order(y)]
    ## test <- sum(resx*resy,na.rm=TRUE)
    ## ahva <- sum(modelx$cumhaz[order(x)]*modely$cumhaz[order(y)],na.rm=TRUE)
    ## test <- test/sqrt(ahva)
    ## pvalues[i] <- 2*pnorm(-abs(test))
    pvalues[i] <- as.numeric(independenceTest(x,y,xs,ys,weight="dabrowska")[3])
    print(mean(pvalues[1:i]))
}

unitest <- sum(qnorm(pvalues))/sqrt(nsim)
2*pnorm(-abs(unitest))

x <- rexp(n)
y <- rexp(n)
modelx <- survfit(Surv(x,xs)~1,ctype=1)
resx <- xs - modelx$cumhaz[order(x)]
modely <- survfit(Surv(y,ys)~1,ctype=1)
resy <- ys - modely$cumhaz[order(y)]
sum(resx*resy)

start <- Sys.time()
result <- independenceTest(x,y,xs,ys)
end <- Sys.time()
difftime(end,start)
result
