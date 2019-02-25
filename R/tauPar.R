tauPar <- function(par = 0, dist = "gamma", output = "tau", type = "alpha"){
    if(dist == "gamma"){
        if(output == "tau"){
            return(par / (2 + par))
        }
        else if(output == "par"){
            return(2*par / (1-par))
        }
        else{
            stop("output has to be either 'tau' or 'par'")
        }
    }
    else if(dist == "posstab"){
        if(type == "alpha"){
            return(1 - par)
        }
        else if(type == "theta"){
            if(output == "tau"){
                return(1 - 1 / par)
            }
            else if(output == "par"){
                return(1 / (1-par))
            }
            else{
                stop("type has to be either 'alpha' or 'theta'")
            }
        }
        else{
            stop("output has to be either 'tau' or 'par'")
        }
    }
    else if(dist == "invgauss"){
        laplace <- function(alpha){
            L <- function(s) exp(alpha - sqrt(alpha) * sqrt(2*s + alpha))
            LL <- function(s){
                exp(alpha) * (-sqrt(alpha)) * (exp(-sqrt(alpha)*sqrt(2*s+alpha))*(-sqrt(alpha)) / (2*s+alpha) -
                                               exp(-sqrt(alpha)*sqrt(2*s+alpha))*(2*s+alpha)^(-3/2))
            }
            f <- function(s) s*L(s)*LL(s)
            4 * integrate(f, 0, Inf)$value - 1
        }
        if(output == "tau"){
            if(type == "alpha"){
                return(laplace(par))
            }
            else if(type == "theta"){
                return(laplace(1/par))
            }
            else{
                stop("type has to be either 'alpha' or 'theta'")
            }
        }
        else if(output == "par"){
            if(type == "alpha"){
                return(uniroot(function(alpha) laplace(alpha) - par, interval = c(.0001, 100))$root)
            }
            else if(type == "theta"){
                return(1 / uniroot(function(alpha) laplace(alpha) - par, interval = c(.0001, 100))$root)
            }
            else{
                stop("type has to be either 'alpha' or 'theta'")
            }
        }
        else{
            stop("output has to be either 'tau' or 'par'")
        }
    }
    else{"dist has to be either 'gamma', 'posstab' or 'invgauss'"}
}
