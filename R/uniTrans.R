#' Transform survival data into uniform marginals using a Cox proportional hazards model
#'
#' @title Transform survival data into uniform marginals using a Cox proportional hazards model
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame in which to interpret the variables named in the
#'          \code{formula} argument.
#' @return Data.frame with two vectors of failure times and two vectors of status indicators
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
uniTrans <- function(formula, data){
    mf <- match.call()
    m <- match(c("formula", "data"), names(mf), 0L)
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if(length(grep("cluster",names(mf))) == 0)
        stop("RHS of formula has to have a 'cluster' term")
    if(ncol(mf)>2){
        suppe <- mf[,2:ncol(mf)]
        suppe <- suppe[,-grep("cluster",names(suppe))]
        m <- coxph(mf[,1]~suppe)
        time <- 1-exp(-predict(m,type="expected"))
    }
    else time <- mf[,1][,1]
    status <- mf[,1][,2]
    clusters <- mf[,grep("cluster",names(mf))]
    firsts <- match(unique(clusters), clusters)
    seconds <- length(clusters) - match(unique(clusters),rev(clusters)) + 1
    if(length(intersect(firsts,seconds))>0){
        firsts <- firsts[!(firsts %in% seconds)]
        seconds <- seconds[!(seconds %in% firsts)]
    }
    if(any(table(clusters)>2))
        warnings("Some clusters have more than 2 observations; these will not be used")
    x <- time[firsts];xstatus <- status[firsts]
    y <- time[seconds];ystatus <- status[seconds]
    data.frame(x,y,xstatus,ystatus)
}
