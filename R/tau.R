tauCens <- function(x,y,xstatus,ystatus, alpha = .05){
  KMx <- prodlim::prodlim(Hist(x, xstatus) ~ 1)
  KMy <- prodlim::prodlim(Hist(y, ystatus) ~ 1)
  tauu <- taucpp(x,y,xstatus,ystatus,KMx$surv,KMy$surv,KMx$time,KMy$time)
  n <- length(x)
  a <- tauu$a
  b <- tauu$b
  var.hat <- 4 * (sum(apply(a, 2, sum)^2) - sum(a^2)) * 
    (sum(apply(b, 2, sum)^2) - sum(b^2)) / (n * (n - 1) * (n - 2)) + 
    2 * sum(a^2) * sum(b^2) / (n * (n - 1))
  var.hat <- var.hat / (sum(a^2) * sum(b^2))
  data.frame(est = tauu$tau, SE = sqrt(var.hat), 
             lwr = tauu$tau-qnorm(1-alpha/2)*sqrt(var.hat), 
             upr = tauu$tau+qnorm(1-alpha/2)*sqrt(var.hat))
}



