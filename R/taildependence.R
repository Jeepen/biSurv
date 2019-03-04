tailDependence <- function(x, y, xstatus, ystatus, q, tail = "lwr", method = "dabrowska"){
  if(method == "dabrowska"){
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    haz <- hazards(x,y,xstatus,ystatus)
    H <- (haz$lambda10 * haz$lambda01 - haz$lambda11) / ((1 - haz$lambda10) * (1 - haz$lambda01))
    H[is.nan(H)] <- 0
    KMx <- prodlim(Hist(x,xstatus) ~ 1)
    KMy <- prodlim(Hist(y,ystatus) ~ 1)
    Fx <- 1 - KMx$surv
    Fy <- 1 - KMy$surv
    if(tail == "lwr"){
      qx <- min(KMx$time[Fx >= q])
      qy <- min(KMy$time[Fy >= q])
      prob <- KMx$surv[KMx$time == qx] * KMy$surv[KMy$time == qy] * prod(1 - H[xuni <= qx, yuni <= qy])
      (2*q - 1 + prob) / q
    }
    else if(tail == "upr"){
      qx <- min(KMx$time[Fx >= q])
      qy <- min(KMy$time[Fy >= q])
      (KMx$surv[KMx$time == qx] * KMy$surv[KMy$time == qy] * prod(1 - H[xuni <= qx, yuni <= qy])) / (1 - q)
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
  }
  else if(method == "fast"){
    KMy <- prodlim(Hist(y, ystatus) ~ 1)
    FF <- 1-KMy$surv
    qq <- min(KMy$time[FF >= q])
    if(tail == "lwr"){
      xx <- x[y < qq & ystatus == 1]
      xxstatus <- xstatus[y < qq & ystatus == 1]
    }
    else if(tail == "upr"){
      xx <- x[y > qq]
      xxstatus <- xstatus[y > qq]
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
    KMx <- prodlim(Hist(xx, xxstatus) ~ 1)
    if(tail == "lwr"){
      Fx <- 1-KMx$surv
      max(Fx[KMx$time < qq])
    }
    else if(tail == "upr"){
      max(KMx$surv[KMx$time > qq])
    }
  }
}

tailDependence <- function(x, y, xstatus, ystatus, q, tail = "lwr", method = "dabrowska", CI = FALSE,
                           nsim = 1000, alpha = .05){
  if(CI == FALSE){
    tailDependence(x, y, xstatus, ystatus, q, tail = tail, method = method)
  }
  else if(CI == TRUE){
    helpfunc <- function(x, y, xstatus, ystatus, q, tail = tail, method = method){
      ord <- sample(1:length(x), replace = T)
      tailDependece(x[ord], y[ord], xstatus[ord], ystatus[ord], q, tail = tail, method = method)
      }
    est <- replicate(nsim, helpfunc(x, y, xstatus, ystatus, q, tail = tail, method = method))
    paste0(tailDependece(x, y, xstatus, ystatus, q, tail = tail, method = method),
             " [", quantile(est, probs = alpha/2), ",", quantile(est, probs = 1-alpha/2),
             "]")
  }
  else{
    stop("CI has to be either 'TRUE' or 'FALSE'")
  }
}