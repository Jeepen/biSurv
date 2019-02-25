CHRres <- function(x,y,xstatus,ystatus, weight = "auto"){
  h <- hazards(x,y,xstatus,ystatus)
  if(weight == "auto"){
    weight <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
  }
  else if(dim(weight) != dim(h$lambda11)){
    stop("dim(weight) has to equal dim()")
  }
  CHRresiduals <- weight * (h$lambda11 - h$lambda10 * h$lambda01)
  CHRresiduals[is.nan(CHRresiduals)] <- 0
  testStatistic <- sum(CHRresiduals)
#  varianceEst <- 
}
