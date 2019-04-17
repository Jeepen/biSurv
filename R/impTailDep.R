#' Theoretical 'tail dependence' at given quantile for given frailty model with given parameter
#'
#' @title Theoretical 'tail dependence' at given quantile for given frailty model with given parameter
#' @param q Quantile
#' @param par Parameter value for frailty model
#' @param dist Frailty distribution. Gamma ("gamma"), positive stable ("posstab") and inverse Gaussian ("invgauss") are available.
#' @param type Parameterization of frailty parameter (either "alpha" or "theta")
#' @param Tail Lower ("lwr") or upper ("upr") tail dependence
#' @return 'Tail dependence' at given quantile for given frailty model with given parameter
#' @seealso tailDep
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
impTailDep <- function(q, par, dist = "gamma", type = "alpha", tail = "lwr"){
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
      (exp(par - (par^2 + 2 * log(1-q) * (log(1-q) - 2 * par))^.5) + 2 * q - 1) / q
    }
    else if(tail == "upr"){
      exp(par - (par^2 + 2*log(1-q) * (log(1-q) - 2*par))^.5) / (1-q)
    }
    else{
      stop("tail has to be either 'lwr' or 'upr'")
    }
  }
  else{
      stop("dist has to be either 'gamma', 'posstab' or 'invgauss'")
  }
}


