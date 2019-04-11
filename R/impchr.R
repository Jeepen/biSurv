CHRtheo <- function(par, dist = "gamma", type = "alpha"){
  if(dist == "gamma"){
    rep(1+par, 97)
  }
  else if(dist == "posstab"){
    if(type == "theta"){par = 1/par}
    1 - (1-alpha) / (alpha * log(seq(.01,.97,.01)))
  }
  else if(dist == "invgauss"){
    if(type == "theta"){par = 1/par}
    1 + 1 / (alpha - log(seq(.01,.97,.01)))
  }
  else{
    stop("dist has to be either 'gamma', 'posstab' or 'invgauss'")
  }
}