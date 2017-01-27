#' @title \emph{tsrToDF()}
#'
#' @description \emph{tsrToDF} is a utility function in \bold{TSRchitect}
#' that converts a \emph{tsrCluster()}-generated list of TSR data to
#' a data frame.
#'
#' @param x A list (of lists) with TSR data as returned by tsrCluster()
#' @return A dataframe of TSRs with variables\cr
#' \enumerate{
#'         \item seq = sequence identifier (seq)
#'         \item start = start of TSR (num)
#'         \item end = end of TSR (num)
#'         \item strand = + or - (factor)
#'         \item nTSSs = count of TSSs (num)
#'         \item tsrWidth = width of TSR (num)
#'         \item shapeIndex = shape index value of TSR (num)
#' }
#' @export


tsrToDF <- function(x) {
    len.list <- length(x)
    seq.vec <- names(x)
    final.matrix <- matrix(NA, nrow=1, ncol=7)

    for (j in 1:len.list) {
        vector(mode="list", length=2) -> my.list.p
        vector(mode="list", length=2) -> my.list.m
        as.character(seq.vec[j]) -> this.seq
        x[[j]] -> this.seq.list

        this.seq.list$plus -> this.p
        c("+") -> my.strand
        matrix(NA, nrow=length(this.p), ncol=7) -> plus.matrix
        if (length(this.p) > 0 ) { # ensures there are + strand TSRs to process
           for (i in 1:max(1,length(this.p))) {
               this.p[[i]] -> my.tsr
               countsToVector(my.tsr) -> my.vec
               shapeIndex(my.vec) -> my.SI
               round(my.SI, digits=2) -> my.SI
               tsrCounts(my.tsr) -> my.counts
               tsrWidth(my.tsr) -> my.width
               range(my.vec) -> my.range
               my.range[1] -> my.start
               my.range[2] -> my.end
               c(this.seq, my.start, my.end, my.strand, my.counts,
                 my.width, my.SI) -> my.string
               my.string -> plus.matrix[i,]
           }
        }

        this.seq.list$minus -> this.m
        c("-") -> my.strand
        matrix(NA, nrow=length(this.m), ncol=7) -> minus.matrix
        if (length(this.m) > 0 ) { # ensures there are - strand TSRs to process
           for (i in 1:max(1,length(this.m))){
               this.m[[i]] -> my.tsr
               countsToVector(my.tsr) -> my.vec
               shapeIndex(my.vec) -> my.SI
               round(my.SI, digits=2) -> my.SI
               tsrCounts(my.tsr) -> my.counts
               tsrWidth(my.tsr) -> my.width
               range(my.vec) -> my.range
               my.range[1] -> my.start
               my.range[2] -> my.end
               c(this.seq, my.start, my.end, my.strand, my.counts,
                 my.width, my.SI) -> my.string
               my.string -> minus.matrix[i,]
           }
        }
 
        rbind(plus.matrix, minus.matrix) -> seq.matrix
        rbind(final.matrix, seq.matrix) -> final.matrix
    }
    final.matrix <- final.matrix[-1,] #removes the 1st row used to initialize
    colnames(final.matrix) <- c("seq", "start", "end", "strand",
                                "nTSSs", "tsrWidth", "shapeIndex")
    final.df <- as.data.frame(final.matrix)
#Convert dataframe column classes to appropriate types:    
    final.df$seq   <- as.character(final.df$seq)
    final.df$start <- as.numeric(as.character(final.df$start))
    final.df$end   <- as.numeric(as.character(final.df$end))
    final.df$nTSSs <- as.numeric(as.character(final.df$nTSSs))
    final.df$tsrWidth   <- as.numeric(as.character(final.df$tsrWidth))
    final.df$shapeIndex   <- as.numeric(as.character(final.df$shapeIndex))
    return(final.df)
}
