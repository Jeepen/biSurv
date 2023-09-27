uniTrans <- function(x,y,,xstatus,ystatus,method="km",X=1,Intercept = TRUE){
    if(method == "km"){
        time<-c(x,y)
        status<-c(xstatus,ystatus)
        km <- prodlim(Hist(time, status) ~ 1)
        predict(km, times = time, type = "cuminc")
    }
    else if(method = "cox"){
        coxph(Surv(time,status)~X)
    }
}








