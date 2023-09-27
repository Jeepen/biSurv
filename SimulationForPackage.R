## Seed&packages&workingDirectory
set.seed(20192102)
setwd('~/Dropbox/RPackage/biSurv')
library(biSurv)
library(foreach)
library(doRNG)
library(doMC)
registerDoMC(4)
library(xtable)
source("~/Dropbox/RPackage/tests/simfunction.R")

## changed function
tauCens0 <- function (formula, data = NULL, method = "adjusted") 
{
    d <- uniTrans(formula, data)
    if (ncol(d) != 4) 
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x
    y <- d$y
    xstatus <- d$xstatus
    ystatus <- d$ystatus
    id <- rep(1:length(x), 2)
    if (method == "adjusted") {
        KMx <- list(time = Rfast::sort_unique(x), surv = biSurv:::KaplanMeier2(x, 
            xstatus))
        KMy <- list(time = Rfast::sort_unique(y), surv = biSurv:::KaplanMeier2(y, 
            ystatus))
        tauu <- biSurv:::taucpp(x, y, xstatus, ystatus, KMx$surv, KMy$surv, 
            KMx$time, KMy$time)
        n <- length(x)
        a <- tauu$a
        b <- tauu$b
        var.hat <- 4 * (sum(apply(a, 2, sum)^2) - sum(a^2)) * 
            (sum(apply(b, 2, sum)^2) - sum(b^2))/(n * (n - 1) * 
            (n - 2)) + 2 * sum(a^2) * sum(b^2)/(n * (n - 1))
        var.hat <- var.hat/(sum(a^2) * sum(b^2))
        tau <- tauu$tau
    }
    else if (method == "naive") {
        obs <- (xstatus == 1 & ystatus == 1)
        xx <- x[obs]
        yy <- y[obs]
        tau <- cor(xx, yy, method = "kendall")
        n <- length(xx)
        var.hat <- (2 * (2 * n + 5))/(9 * n * (n - 1))
    }
    else {
        stop("method has to be either 'adjusted' or 'naive'")
    }
    se <- sqrt(var.hat)
    out <- list(tau = tau, se = se)
    out
}

## Parameters and initialization
taus <- c(0,.25,.5, .75)
sims <- 10000
results <- list(clayton = matrix(NA, nrow = sims, ncol = length(taus)), gumbel = matrix(NA, nrow = sims, ncol = length(taus)),
                claytonSE = matrix(NA, nrow = sims, ncol = length(taus)), gumbelSE = matrix(NA, nrow = sims, ncol = length(taus)))
colnames(results$clayton) <- colnames(results$gumbel) <- colnames(results$claytonSE) <- colnames(results$gumbelSE) <- taus

## Simulation
startTime <- Sys.time()
results <- foreach(i = 1:sims, .options.RNG = 20022021) %dorng% {
    if(i %% 100 == 0)
        print(i)
    clayton <- gumbel <- list()
    for(j in 1:length(taus)){
        simClayton <- simFrail(100, par = tauPar(taus[j], output = "par"))
        timecens <- rexp(200, .5)
        simClayton <- transform(simClayton, time = pmin(time, timecens), status = as.numeric(timecens > time))
        simGumbel <- simFrail(100, dist = "posstab", par = tauPar(taus[j], dist = "posstab", output = "par", type = "alpha"))
        simGumbel <- transform(simGumbel, time = pmin(time, timecens), status = as.numeric(timecens > time))
        clayton[[j]] <- tauCens0(Surv(time,status) ~ cluster(id), data = simClayton)
        gumbel[[j]] <- tauCens0(Surv(time, status) ~ cluster(id), data = simGumbel)
    }
    if(i %% 100 == 0)
        print(list(clayton = clayton, gumbel = gumbel))
    list(clayton = clayton, gumbel = gumbel)
}
difftime(Sys.time(), startTime)

table <- rbind(
    sapply(1:4, function(i) mean(do.call("c", lapply(results, function(x) x$clayton[[i]]$tau)))),
    sapply(1:4, function(i) mean(do.call("c", lapply(results, function(x) x$clayton[[i]]$se)))),
    sapply(1:4, function(i) sd(do.call("c", lapply(results, function(x) x$clayton[[i]]$tau)))),
    sapply(1:4, function(i) mean(do.call("c", lapply(results, function(x) x$gumbel[[i]]$tau)))),
    sapply(1:4, function(i) mean(do.call("c", lapply(results, function(x) x$gumbel[[i]]$se)))),
    sapply(1:4, function(i) sd(do.call("c", lapply(results, function(x) x$gumbel[[i]]$tau))))
)
rownames(table) <- c("Clayton", "Clayton theoretical SE", "Clayton empirical SE", "Gumbel", "Gumbel theoretical SE","Gumbel empirical SE")
xtable(table)

