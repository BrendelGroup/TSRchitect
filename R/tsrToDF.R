# An internal function for TSS clustering
# Creates a data frame for a given output of .tsrCluster

tsrToDF <- function(x) {
    len.list <- length(x)
    final.matrix <- matrix(NA, nrow=1, ncol=4)
    for (j in 1:len.list) {
        vector(mode="list", length=2) -> my.list.p
        vector(mode="list", length=2) -> my.list.m
        as.character(chr.vec[j]) -> this.chr
        x[[j]] -> this.chr.list
        this.chr.list$plus -> this.p
        c("+") -> my.strand
        matrix(NA, nrow=length(this.p), ncol=4) -> plus.matrix

        for (i in 1:length(this.p)) {
            this.p[[i]] -> my.vec
            range(my.vec) -> my.range
            my.range[1] -> my.start
            my.range[2] -> my.end
            c(this.chr, my.start, my.end, my.strand) -> my.string
            my.string -> plus.matrix[i,]
        }
            this.chr.list$minus -> this.m
            c("-") -> my.strand
            matrix(NA, nrow=length(this.m), ncol=4) -> minus.matrix

        for (i in 1:length(this.m)){
            this.m[[i]] -> my.vec
            range(my.vec) -> my.range
            my.range[1] -> my.start
            my.range[2] -> my.end
            c(this.chr, my.start, my.end, my.strand) -> my.string
            my.string -> minus.matrix[i,]
        }
        rbind(plus.matrix, minus.matrix) -> chr.matrix
    }
    final.matrix <- rbind(final.matrix, chr.matrix)
    colnames(final.matrix) <- c("chr", "start", "end", "strand")
    final.df <- as.data.frame(final.matrix)
    return(final.df)
}
    
    
