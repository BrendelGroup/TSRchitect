#' An internal function for TSS clustering
#' Creates a data frame for a given output of .tsrCluster
#' @export


tsrToDF <- function(x) {
    len.list <- length(x)
    chr.vec <- names(x)
    final.matrix <- matrix(NA, nrow=1, ncol=7)

    for (j in 1:len.list) {
        vector(mode="list", length=2) -> my.list.p
        vector(mode="list", length=2) -> my.list.m
        as.character(chr.vec[j]) -> this.chr
        x[[j]] -> this.chr.list
        this.chr.list$plus -> this.p
        c("+") -> my.strand
        matrix(NA, nrow=length(this.p), ncol=7) -> plus.matrix

        for (i in 1:length(this.p)) {
            this.p[[i]] -> my.tsr
            countsToVector(my.tsr) -> my.vec
            shapeIndex(my.vec) -> my.SI
            round(my.SI, digits=2) -> my.SI
            tsrCounts(my.tsr) -> my.counts
            tsrWidth(my.tsr) -> my.width
            range(my.vec) -> my.range
            my.range[1] -> my.start
            my.range[2]+1 -> my.end
            c(this.chr, my.start, my.end, my.strand, my.counts, my.width, my.SI) -> my.string
            my.string -> plus.matrix[i,]
        }
            this.chr.list$minus -> this.m
            c("-") -> my.strand
            matrix(NA, nrow=length(this.m), ncol=7) -> minus.matrix

        for (i in 1:length(this.m)){
            this.m[[i]] -> my.tsr
            countsToVector(my.tsr) -> my.vec
            shapeIndex(my.vec) -> my.SI
            round(my.SI, digits=2) -> my.SI
            tsrCounts(my.tsr) -> my.counts
            tsrWidth(my.tsr) -> my.width
            range(my.vec) -> my.range
            my.range[1]-1 -> my.start
            my.range[2] -> my.end
            c(this.chr, my.start, my.end, my.strand, my.counts, my.width, my.SI) -> my.string
            my.string -> minus.matrix[i,]
        }
    rbind(plus.matrix, minus.matrix) -> chr.matrix
    rbind(final.matrix, chr.matrix) -> final.matrix
    }
    final.matrix <- final.matrix[-1,] #removes the empty first row used to initilize the matrix
    colnames(final.matrix) <- c("chr", "start", "end", "strand", "nTSSs", "tsrWidth", "shapeIndex")
    final.df <- as.data.frame(final.matrix)
    return(final.df)
}


