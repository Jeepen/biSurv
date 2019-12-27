#' Transform survival data into uniform marginals using a Cox proportional hazards model
#'
#' @title Transform survival data into uniform marginals using a Cox proportional hazards model
#' @param formula a formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term
#' @param data a data.frame containing the variables in the model
#' @param cluster Cluster variable
#' @return Data.frame with two vectors of failure times and two vectors of status indicators
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
uniTrans <- function(formula, data = NULL, cluster = NULL){
    mf <- match.call()
    cluster <- eval(substitute(cluster),data)
    m <- match(c("formula", "data", "cluster"), names(mf), 0L)
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if(ncol(model.matrix(formula, data))>1){
        m <- coxph(formula, data)
        time <- 1-exp(-predict(m,type="expected"))
    }
    else time <- mf[,1][,1]
    status <- mf[,1][,2]
    if(is.null(cluster)){
        data.frame(time, status)
    }
    else{
        firsts <- match(unique(cluster), cluster)
        seconds <- length(cluster) - match(unique(cluster),rev(cluster)) + 1
        if(length(intersect(firsts,seconds))>0){
            firsts <- firsts[!(firsts %in% seconds)]
            seconds <- seconds[!(seconds %in% firsts)]
        }
        if(any(table(cluster)>2))
            warning("Some clusters have more than 2 observations; these will not be used")
        x <- time[firsts];xstatus <- status[firsts]
        y <- time[seconds];ystatus <- status[seconds]
        data.frame(x,y,xstatus,ystatus)
    }
}
