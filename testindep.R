library(biSurv)

nsim <- 1000
pvalues <- numeric(nsim)
for(i in 1:nsim){
  if(i%%10 == 0){print(i)}
  xstar <- rexp(100)
  ystar <- rexp(100)
  xc <- rexp(100, .3)
  yc <- rexp(100, .3)
  x <- pmin(xstar, xc)
  y <- pmin(ystar, yc)
  xstatus <- (xstar <= xc)
  ystatus <- (ystar <= yc)
  pvalues[i] <- CHRres(x,y,xstatus,ystatus)$`p-value`
}
summary(pvalues)
mean(pvalues < 0.05)
hist(pvalues)
