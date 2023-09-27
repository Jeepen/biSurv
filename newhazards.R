hazard2 <- function(x,y,xstatus,ystatus, xuni = sort_unique(x), yuni = sort_unique(y)){
    R <- do.call("rbind", lapply(xuni, FUN = function(xuni){ 
        sapply(yuni, FUN = function(yuni) sum(x >= xuni & y >=yuni))}))
    L11 <- do.call("rbind", lapply(xuni, FUN = function(xuni){
        sapply(yuni, FUN = function(yuni){sum(xstatus*ystatus*(x==xuni&y==yuni))})})) / R
    L10 <- do.call("rbind", lapply(xuni, FUN = function(xuni){
        sapply(yuni, FUN = function(yuni){sum(xstatus*(x==xuni&y>=yuni))})})) / R
    L01 <- do.call("rbind", lapply(xuni, FUN = function(xuni){
        sapply(yuni, FUN = function(yuni){sum(ystatus*(x>=xuni&y==yuni))})})) / R
    list(lambda11 = L11, lambda10 = L10, lambda01 = L01)
}
