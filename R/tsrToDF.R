#' @title \emph{tsrToDF()}
#'
#' @description \emph{tsrToDF} is an internal utility function in
#' that converts a \emph{tsrCluster()}-generated list of TSR data to
#' a data frame.
#'
#' @param x A list with TSR data as returned by tsrCluster()
#'
#' @keywords internal
#'
#' @return A data frame of TSRs with variables\cr
#' \enumerate{
#'         \item seq = sequence identifier (seq)
#'         \item start = start of TSR (num)
#'         \item end = end of TSR (num)
#'         \item strand = + or - (factor)
#'         \item nTSSs = count of TSSs (num)
#'         \item tsrWidth = width of TSR (num)
#'         \item shapeIndex = shape index value of TSR (num)
#' }

tsrToDF <- function(clstOut) {
    
                    df.fun.p <- function(chrLst) {
                          my.string <- lapply(chrLst, countsToVector)
                          my.SI <- sapply(my.string, shapeIndex)
                          my.SI <- round(my.SI, digits=2)
                          my.vec <- lapply(chrLst, "[", 1, )
                          my.range <- lapply(my.vec, range)
                          my.width <- lapply(my.range, function(x) (x[2]-x[1])+1)
                          my.width <- do.call(rbind, my.width)
                          my.start <- vapply(my.range, min, numeric(1))
                          my.end <- vapply(my.range, max, numeric(1))
                          my.counts <- vapply(chrLst, tsrCounts, numeric(1))
                          my.strand <- rep("+", length(my.vec))
                          string.out.p <- cbind(my.start, my.end, my.strand, my.counts, my.width, my.SI)
                          return(string.out.p)
                      }

                     df.fun.m <- function(chrLst) {
                          my.string <- lapply(chrLst, countsToVector)
                          my.SI <- sapply(my.string, shapeIndex)
                          my.SI <- round(my.SI, digits=2)
                          my.vec <- lapply(chrLst, "[", 1, )
                          my.range <- lapply(my.vec, range)
                          my.width <- lapply(my.range, function(x) (x[2]-x[1])+1)
                          my.width <- do.call(rbind, my.width)
                          my.start <- vapply(my.range, min, numeric(1))
                          my.end <- vapply(my.range, max, numeric(1))
                          my.counts <- vapply(chrLst, tsrCounts, numeric(1))
                          my.strand <- rep("-", length(my.vec))
                          string.out.m <- cbind(my.start, my.end, my.strand, my.counts, my.width, my.SI)
                          return(string.out.m)
                      }

                    chrToDF  <- function(x) {
                        my.list <- vector(mode="list")
                        list.p <- x$plus
                        list.m <- x$minus
                        my.list$plus <- df.fun.p(list.p)
                        my.list$minus <- df.fun.m(list.m)
                        my.out <- do.call(rbind, my.list)
                        return(my.out)
                    }
                    
    options(warn=-1)
    names.list <- lapply(clstOut, function(x) sapply(x, length))
    my.sum <- lapply(names.list, sum)
    seq <- as.data.frame(rep(names(my.sum), my.sum))
    this.list <- lapply(clstOut, chrToDF)
    my.df <- do.call(rbind, this.list)
    my.df <- as.data.frame(my.df)
    my.df <- data.frame(seq, my.df)
    colnames(my.df) <- c("seq", "start","end", "strand", "nTSSs", "tsrWidth", "shapeIndex")
    return(my.df)
}


