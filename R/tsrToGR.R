# An internal function for TSS clustering
# Creates a GRanges object for a given output of tsrCluster

tsrToGR <- function(x, seqname="chr1") {
    my.list <- vector(mode="list",length=2)
    my.list$plus <- vector(mode="list")
    my.list$minus <- vector(mode="list")
    this.p <-  x$plus
    my.strand <- c("+")
    for (i in 1:length(this.p)) {
        this.p[[i]] -> my.vec
        range(my.vec) -> my.range
        my.range[1] -> my.start
        my.range[2] -> my.end
        GRanges(seqnames = seqname, strand = c(my.strand),
                 ranges = IRanges(start = c(my.start), end=my.end)
                ) -> myGR
        myGR -> my.list$plus[[i]]
    }
    this.m <- x$minus
    my.strand <- c("-")
    for (i in 1:length(this.m)){
        this.m[[i]] -> my.vec
        range(my.vec) -> my.range
        my.range[1] -> my.start
        my.range[2] -> my.end
        GRanges(seqnames = seqname, strand = c(my.strand),
                 ranges = IRanges(start = c(my.start), end=my.end)
                ) -> myGR
        myGR -> my.list$minus[[i]]
    }
    return(my.list)
}
    
    
