CHRres <- function(x,y,xstatus,ystatus, weight = "auto"){
  h <- hazards(x,y,xstatus,ystatus)
  if(weight == "auto"){
    weight <- eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
  }
  else if(dim(weight) != dim(h$lambda11)){
    stop("dim(weight) has to equal dim()")
  }
  else{
    eyy <- eyyfunc(x,y,sort_unique(x),sort_unique(y))
  }
  CHRresiduals <- weight * (h$lambda11 - h$lambda10 * h$lambda01)
  CHRresiduals[is.nan(CHRresiduals)] <- 0
  testStatistic <- -abs(sum(CHRresiduals))
  varianceEst <- sum(weight^2 * eyy * h$lambda10 * h$lambda01)
  test <- testStatistic / sqrt(varianceEst)
  d <- data.frame(testStatistic, test, 2*pnorm(test))
  colnames(d) <- c("Test statistic", "Standardized test statistic", "p value")
  d
}
