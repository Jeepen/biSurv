EMsurv2 <- function(formula, data, tol = .001, maxIt = 10){
    Call <- match.call()
    d <- uniTrans(formula, data)
    if(ncol(d) != 4) 
        stop("RHS needs a 'cluster(id)' element")
    x <- d$x; y <- d$y; xstatus <- d$xstatus; ystatus <- d$ystatus
    xuni <- sort_unique(x)
    yuni <- sort_unique(y)
    tmp <- EMcpp(x,y,xstatus,ystatus,xuni,yuni, tol, maxIt)
    n1 <- length(xuni); n2 <- length(yuni)
    t(apply(apply(tmp[n1:1,n2:1], 2, cumsum), 
    1, cumsum))[n1:1,n2:1]
}
