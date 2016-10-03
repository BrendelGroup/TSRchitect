# An internal function for TSS clustering
# Creates a GRanges object for a given output of tsrCluster

tsrToGR <- function(x) {
    len.list <- length(x)
    chr.vec <- names(x)
    all.list <- vector(mode="list")
    my.list <- vector(mode="list")
    for (j in 1:len.list) {
        vector(mode="list", length=2) -> my.list.p
        vector(mode="list", length=2) -> my.list.m
        as.character(chr.vec[j]) -> this.chr
        x[[j]] -> this.chr.list
        this.chr.list$plus -> this.p
        my.strand <- c("+")
        for (i in 1:length(this.p)) {
            this.p[[i]] -> my.vec
            range(my.vec) -> my.range
            my.range[1] -> my.start
            my.range[2] -> my.end
            GRanges(seqnames = this.chr, strand = c(my.strand),
                    ranges = IRanges(start = c(my.start), end=my.end)
                    ) -> myGR
            is(myGR)
            myGR -> my.list.p[[i]]
        }
            this.chr.list$minus -> this.m
            my.strand <- c("-")
        for (i in 1:length(this.m)){
            this.m[[i]] -> my.vec
            range(my.vec) -> my.range
            my.range[1] -> my.start
            my.range[2] -> my.end
            GRanges(seqnames = this.chr, strand = c(my.strand),
                 ranges = IRanges(start = c(my.start), end=my.end)
                ) -> myGR
            myGR -> my.list.m[[i]]
        }
        my.list.p -> my.list$plus
        my.list.m -> my.list$minus
        my.list -> all.list[[j]]
    }
    
    names(all.list) <- chr.vec
    return(my.list)
}
    
    
