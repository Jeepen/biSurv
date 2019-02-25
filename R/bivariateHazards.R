hazards <- function(x,y,xstatus,ystatus){
  xuni <- sort_unique(x)
  yuni <- sort_unique(y)
  hazardscpp(x,y,xstatus,ystatus,xuni,yuni)
}
