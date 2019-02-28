impliedTailDependence <- function(q, par, dist = "gamma", type = "alpha", tail = "lwr"){
  if(dist == "gamma"){
    if(tail == "lwr"){
      ((2*(1-q)^(-par) - 1)^(-1/par) + 2*q - 1) / q
    }
    else if(tail == "upr"){
      (2*(1-q)^(-par)-1)^(-1/par)/ (1-q)
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
  }
  else if(dist == "posstab"){
    if(type == "theta"){
      par <- 1/par
    }
    if(tail == "lwr"){
      (exp(-(2*(-log(1-q))^(1/par))^par) + 2*q - 1) / q
    }
    else if(tail == "upr"){
      exp(-(2*(-log(1-q))^(1/par))^par)/ (1-q)
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
  }
  else if(dist == "invgauss"){
    if(type == "theta"){
      par <- 1/par
    }
    if(tail == "lwr"){
      invgausstheolwrMZ <- function(q) (exp(par - (par^2 + 2 * log(1-q) * (log(1-q) - 2 * par))^.5) + 2 * q - 1) / q
    }
    else if(tail == "upr"){
      exp(par - (par^2 + 2*log(1-q) * (log(1-q) - 2*par))^.5) / (1-q)
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
  }
}


