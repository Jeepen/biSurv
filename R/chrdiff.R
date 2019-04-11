chrdiff <- function(x, y, xstatus, ystatus, par, dist = "gamma", type = "alpha", n = 5){
  imp <- CHR(x,y,xstatus,ystatus, n = n)
  imp0 <- numeric(97)
  breaks <- seq(0,1,length.out=n+1)
  for(i in 1:(n-1)){
    imp0[(1+(i-1)*100/n):(i*100/n)] <- imp[i]
  }
  imp0[(1+(n-1)*100/n):97] <- imp[n]
  theo <- CHRtheo(par, dist = disk, type = type)
  sum((imp0 - theo)^2)
}
