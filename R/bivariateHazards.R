hazards <- function(x,y,xstatus,ystatus){
  xuni <- sort_unique(x)
  yuni <- sort_unique(y)
  h <- hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
  h$lambda11[is.nan(h$lambda11)] <- 0
  h$lambda10[is.nan(h$lambda10)] <- 0
  h$lambda01[is.nan(h$lambda01)] <- 0
  h
}
