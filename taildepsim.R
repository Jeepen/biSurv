## Seed&packages&workingDirectory
set.seed(20192102)
setwd('/home/jeep/Dropbox/RPackage/biSurv')
library(biSurv)
library(copula)
library(xtable)

## Parameters and initialization
taus <- c(0,.25,.5, .75)
sims <- 10000
results <- list(lwr = array(NA, dim = c(sims, length(taus), 4), dimnames = list(1:sims, taus, c("clayton", "claytonFast", "Gumbel", "GumbelFast"))),
                upr = array(NA, dim = c(sims, length(taus), 4), dimnames = list(1:sims, taus, c("clayton", "claytonFast", "Gumbel", "GumbelFast"))))


## Simulation
startTime <- Sys.time()
for(i in 1:sims){
    if(i %% 100 == 0){
        print(i)
    }
    for(j in 1:length(taus)){
        cop.clayton <- claytonCopula(param = tauPar(taus[j], output = "par"), dim = 2, use.indepC = "TRUE")
        simClayton <- rCopula(100, cop.clayton)
        xstar <- qexp(simClayton[,1])
        ystar <- qexp(simClayton[,2])
        xc <- rexp(100, .3)
        yc <- rexp(100, .3)
        x <- pmin(xstar, xc)
        y <- pmin(ystar, yc)
        xstatus <- (xstar <= xc)
        ystatus <- (ystar <= yc)
        cop.gumbel <- gumbelCopula(tauPar(taus[j], dist = "posstab", output = "par", type = "theta"), dim = 2, use.indepC="TRUE")
        simGumbel <- rCopula(100, cop.gumbel)
        xstarg <- qexp(simGumbel[,1])
        ystarg <- qexp(simGumbel[,2])
        xcg <- rexp(100, .3)
        ycg <- rexp(100, .3)
        xg <- pmin(xstarg, xcg)
        yg <- pmin(ystarg, ycg)
        xstatusg <- (xstarg <= xcg)
        ystatusg <- (ystarg <= ycg)
        ## Lower tail
        results$lwr[i,j,1] <- tailDependence(x, y, xstatus, ystatus, q = .1)
        results$lwr[i,j,2] <- tailDependence(x, y, xstatus, ystatus, q = .1, method = "fast")
        results$lwr[i,j,3] <- tailDependence(xg, yg, xstatusg, ystatusg, q = .1)
        results$lwr[i,j,4] <- tailDependence(xg, yg, xstatusg, ystatusg, q = .1, method = "fast")
        ## Upper tail
        results$upr[i,j,1] <- tailDependence(x, y, xstatus, ystatus, q = .8, tail = "upr")
        results$upr[i,j,2] <- tailDependence(x, y, xstatus, ystatus, q = .8, method = "fast", tail = "upr")
        results$upr[i,j,3] <- tailDependence(xg, yg, xstatusg, ystatusg, q = .8, tail = "upr")
        results$upr[i,j,4] <- tailDependence(xg, yg, xstatusg, ystatusg, q = .8, method = "fast", tail = "upr")
    }
}
difftime(Sys.time(), startTime)
results$lwr[results$lwr == -Inf] <- NA
results$upr[results$upr == -Inf] <- NA
xtable(t(apply(results$lwr, c(2,3), mean, na.rm = TRUE)), digits = 3)
xtable(t(apply(results$upr, c(2,3), mean, na.rm = TRUE)), digits = 3)
