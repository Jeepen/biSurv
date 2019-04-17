#' Get theoretical CHR as a function of the survival function for different frailty models
#'
#' @title Get theoretical CHR as a function of the survival function for different frailty models
#' @param par Parameter value for frailty model
#' @param dist Frailty distribution. Gamma ("gamma"), positive stable ("posstab") and inverse Gaussian ("invgauss") are available.
#' @return A vector of length 97 with the theoretical CHR when the survival function is equal to 0.01,...,0.97
#' @seealso chrdiff
#' @export
#' @author Jeppe E. H. Madsen
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
