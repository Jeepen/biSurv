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
    if(tail == "lwr"){
      if(type == "theta"){
        par <- 1/par
      }
      (exp(-(2*(-log(1-q))^(1/par))^par) + 2*q - 1) / q
    }
  }
}



ClaytontheolwrMZ <- function(q) ((2*(1-q)^(-.4774438) - 1)^(-1/.4774438) + 
                                   2*q - 1) / q
GumbeltheolwrMZ <- function(q) (exp(-(2*(-log(1-q))^1.162024)^0.8605675) + 
                                  2*q - 1) / q 
invgausstheolwrMZ <- function(q) (exp(0.9665058 - (0.9665058^2 + 2 * log(1-q) *
                                                     (log(1-q) - 2 * 0.9665058))^.5) + 
                                    2 * q - 1) / q


ClaytontheoMZ <- function(q) (2*(1-q)^(-.4774438)-1)^(-1/.4774438)/ (1-q)
GumbeltheoMZ <- function(q) exp(-(2*(-log(1-q))^1.162024)^0.8605675)/ (1-q)
invgausstheoMZ <- function(q) exp(0.9665058-(0.9665058^2 + 2*log(1-q)*
                                               (log(1-q) -2*0.9665058))^.5) / (1-q)