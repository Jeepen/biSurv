#' Transform survival data into uniform marginals using a Cox proportional hazards model
#'
#' @title Transform survival data into uniform marginals using a Cox proportional hazards model
#' @param formula A formula object, with the response on the left of a ~
#'          operator, and the terms on the right.  The response must be a
#'          survival object as returned by the \code{Surv} function. The RHS must contain a 'cluster' term.
#' @param data A data.frame containing the variables in the model.
#' @details Transforms event times to uniform marginals using a Cox model. This is convenient if
#' you want to analyze your data taking covariates into account. If formula doesn't contain any
#' covariates, the function just returns the original data split up in different columns. 
#' @return Data.frame with two vectors of failure times and two vectors of status indicators
#' @export
#' @author Jeppe E. H. Madsen <jeppe.ekstrand.halkjaer@gmail.com>
uniTrans <- function(formula, data = NULL){
    mf <- match.call()
    Terms <- terms(formula, "cluster", data = data)
    mf$formula <- Terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    m <- grep("cluster", names(mf))
    if(length(m) != 0) cluster <- mf[,m] else cluster <- NULL
    X <- model.matrix(formula, data)
    hm <- grep("cluster", colnames(X))
    X <- X[,-hm]
    Y <- model.extract(mf, "response")
    if (!is.Surv(Y)) 
        stop("Response must be a survival object")
    if(any(X != 1)){
        m <- coxph(Y ~ X)
        time <- 1-exp(-predict(m,type="expected"))
    }
    else time <- Y[,1]
    status <- Y[,2]
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
