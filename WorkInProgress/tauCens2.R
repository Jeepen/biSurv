tauCens2 <- function (formula, data = NULL, alpha = 0.05, method = "adjusted") 
{
    Call <- match.call()
    gamma <- emfrail(formula = formula, data = data)
    theta <- 1/exp(gamma$logtheta)
    stable <- emfrail(formula = formula, data = data, distribution = emfrail_dist(dist = "stable"))
    alpha1 <- exp(stable$logtheta)/(1 + exp(stable$logtheta))
    invgauss <- emfrail(formula = formula, data = data, distribution = emfrail_dist(dist = "pvf"))
    alpha2 <- exp(invgauss$logtheta)
    models <- c(tauPar(theta), tauPar(alpha1, dist = "posstab"), 
        tauPar(alpha2, dist = "invgauss"))
    d <- uniTrans(formula, data)
    if (ncol(d) != 4) 
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x
    y <- d$y
    xstatus <- d$xstatus
    ystatus <- d$ystatus
    id <- rep(1:length(x), 2)
    if (method == "adjusted") {
        KMx <- list(time = sort_unique(x), surv = KaplanMeier(x, 
            xstatus))
        KMy <- list(time = sort_unique(y), surv = KaplanMeier(y, 
            ystatus))
        tauu <- taucpp2(x, y, xstatus, ystatus, KMx$surv, KMy$surv, 
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
    se <- c(sqrt(var.hat), sqrt(gamma$var_logtheta) * 2 * exp(gamma$logtheta)/((1 + 
        2 * exp(gamma$logtheta))^2), sqrt(stable$var_logtheta) * 
        abs(exp(stable$logtheta)/(1 + exp(stable$logtheta)) - 
            exp(2 * stable$logtheta)/((1 + exp(stable$logtheta))^2)), 
        NA)
    out <- list(call = Call, coefficients = c(tau, models), se = se)
    class(out) <- "tauCens"
    out
}

KaplanMeier <- function(time, status){
    t <- sort_unique(time)
    Y <- sapply(t, function(x) sum(time>=x))
    dN <- sapply(t, function(x) sum(time==x & status == 1))
    cumprod(1-dN/Y)
}
