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


tsrToDF <- function(x) {
    len.list <- length(x)
    seq.vec <- names(x)
    final.matrix <- matrix(NA, nrow=1, ncol=7)

    for (j in 1:len.list) { ## replace this loop with do.call(rbind,)
        my.list.p <- vector(mode="list", length=2)
        my.list.m <- vector(mode="list", length=2)
        this.seq <- as.character(seq.vec[j])
        this.seq.list <- x[[j]]

        this.p <- this.seq.list$plus
        my.strand <- c("+")
        plus.matrix <- matrix(NA, nrow=length(this.p), ncol=7) #rather than matrix, comput the columns in vectorized fashion
        if (length(this.p) > 0 ) { # ensures there are + strand TSRs to process
           for (i in 1:max(1,length(this.p))) { #this is probaly vectorizable using vapply() #this data 
               my.tsr <- this.p[[i]] 
               my.vec <- countsToVector(my.tsr)
               my.SI <- shapeIndex(my.vec)
               my.SI <- round(my.SI, digits=2)
               my.counts <- tsrCounts(my.tsr)
               my.width <- tsrWidth(my.tsr)
               my.range <- range(my.vec)
               my.start <- my.range[1]
               my.end <- my.range[2]
               my.string <- c(this.seq, my.start, my.end, my.strand, my.counts,
                              my.width, my.SI)
               plus.matrix[i,] <- my.string
           }
        }

        this.m <- this.seq.list$minus
        my.strand <- c("-")
        minus.matrix <- matrix(NA, nrow=length(this.m), ncol=7)
        if (length(this.m) > 0 ) { # ensures there are - strand TSRs to process
           for (i in 1:max(1,length(this.m))){
               my.tsr <- this.m[[i]]
               my.vec <- countsToVector(my.tsr)
               my.SI <- shapeIndex(my.vec)
               my.SI <- round(my.SI, digits=2)
               my.counts <- tsrCounts(my.tsr)
               my.width <- tsrWidth(my.tsr)
               my.range <- range(my.vec)
               my.start <- my.range[1]
               my.end <- my.range[2]
               my.string <- c(this.seq, my.start, my.end, my.strand, my.counts,
                 my.width, my.SI)
               minus.matrix[i,] <- my.string
           }
        }

        seq.matrix <- rbind(plus.matrix, minus.matrix)
        final.matrix <- rbind(final.matrix, seq.matrix)
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
