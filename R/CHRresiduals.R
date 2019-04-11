CHRres <- function(x,y,xstatus,ystatus, weight = "auto"){
  h <- hazards(x,y,xstatus,ystatus)
  if(weight == "auto"){
    weight <- eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
  }
  else if(dim(weight) != dim(h$lambda11)){
    stop("dim(weight) has to equal dim(lambda11)")
  }
  else{
    eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
  }
  CHRresiduals <- weight * (h$lambda11 - h$lambda10 * h$lambda01)
  CHRresiduals[is.nan(CHRresiduals)] <- 0
  testStatistic <- sum(CHRresiduals)
  varianceEst <- sum(weight^2 * eyy * h$lambda10 * h$lambda01) / length(x)
  test <- -abs(testStatistic) / sqrt(varianceEst * length(x))
  d <- data.frame(testStatistic, sqrt(varianceEst), 2*pnorm(test))
  colnames(d) <- c("Test statistic", "SE", "p-value")
  d
}
