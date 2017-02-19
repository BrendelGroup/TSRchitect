#' @title countsToVector
#' @description An internal function that extracts the values
#' from the output returned from tsrCluster
#'
#' @keywords internal

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

#' @title tsrCounts
#' @description An internal function that caculates the total number of
#' counts associated with a given TSR
#'
#' @keywords internal
#'
#' @return Returns the number of TSS tags associated with a TSR

tsrCounts <- function(x) {
    count.vec <- x[2,]
    my.sum <- sum(count.vec)
    return(my.sum)
}

#' @title tsrWidth
#' @description An internal function that caculates the width
#' of a given TSR from the output of tsrCluster
#'
#' @keywords internal
#'
#' @return Returns a width value (in bp) for a given TSR

tsrWidth <- function(x) {
    coord.vec <- x[1,]
    my.range <- range(coord.vec)
    my.width <- abs(my.range[2]-my.range[1])
    my.width <- my.width+1 #the minimum width possible is 1, not 0
    return(my.width)
}

#' @title shapeIndex
#' @description An internal function that caculates the shape index (SI)
#' of a given TSR from the output of tsrCluster
#'
#' @keywords internal
#'
#' @return Calculates the shape index (SI) for a given TSR

shapeIndex <- function(x) {
        total.size <- length(x)
        unique.tss <- unique(x)
        n.unique.tss <- length(unique.tss)
        p.array <- array(NA,c(1,n.unique.tss))
        for (i in 1:n.unique.tss) {
            n.o <- x
            n.i <- length(which(unique.tss[i]==n.o))
            p.sub.i <- n.i/total.size
            p.log.2 <- log2(p.sub.i)
            p.array[1,i] <- (p.sub.i*p.log.2)
        }
        p.total <- sum(p.array)
        SI <- 2+p.total
        return(SI)
    }
