#' An internal function that extracts the values from the row of counts returned by .tsrCluster
#' @export

countsToVector <- function(x) {
    coord.vec <- x[1,]
    count.vec <- x[2,]
    rep.vec  <- vector(mode="numeric", length=0)
    
    for (i in 1:length(count.vec)) {
        count.vec[i] -> this.count
        coord.vec[i] -> this.coord
        rep(this.coord, this.count) -> this.vec
        c(rep.vec, this.vec) -> rep.vec
        }

    names(rep.vec) <- NULL
    return(rep.vec)
}

#' An internal function that caculates the total number of counts associated with a given TSR
#' @export

tsrCounts <- function(x) {
    count.vec <- x[2,]
    my.sum <- sum(count.vec)
    return(my.sum)
}

#' An internal function that caculates the width of a given TSR
#' @export

tsrWidth <- function(x) {
    coord.vec <- x[1,]
    my.range <- range(coord.vec)
    my.width <- abs(my.range[2]-my.range[1])
    return(my.width)
}
