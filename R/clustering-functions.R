################################################################################
#' tagCountTSS
#' @description an internal function that returns a matrix [a, h] where a = the
#' number of unique TSSs and h = the # of tags observed at that position
#'
#' @import BiocGenerics
#' @import GenomicRanges
#' @importFrom utils write.table
#'
#' @keywords internal
#' @return a matrix [a, h] containing the number of unique TSSs (a) and their
#' abundances (h).


tagCountTSS <- function(y, outfname="TSS.txt", writeDF=FALSE) {
    x <- S4Vectors::split(y,seqnames(y))
    n.seq <- length(x)

    my.matrix <- NULL
    for (i in 1:n.seq) {
#VB Note: Print a progress note on every 20th sequence; 20 should be a parameter
        if (i%%20 == 0) {
            message("... tagCountTSS running with sequence ", i,
                " of ", n.seq, " for TSS set ", outfname, "\n")
        }
        this.seq <- as.character(x[[i]]@seqnames@values[1])

        #starting with the plus strand:
        tss.vec <- start(x[[i]][strand(x[[i]]) == "+"])
        if (length(tss.vec) > 3) {        #stop if there are nearly no tags
            my.TSSs <- unique(tss.vec)
            my.matrix.p <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

            this.TSS <- tss.vec[1]
            n.TSSs <- 1
            k <- 0
            for (j in 2:length(tss.vec)) {
                if (tss.vec[j] == this.TSS) {
                    n.TSSs <- n.TSSs + 1
                }
                else {
                    k <- k + 1
                    my.matrix.p[k,] <- c(this.seq, this.TSS, "+", n.TSSs)
                    this.TSS <- tss.vec[j]
                    n.TSSs <- 1
                }
            }
            k <- k + 1
            my.matrix.p[k,] <- c(this.seq, this.TSS, "+", n.TSSs)
# ... add the plus strand matrix of this.seq to the overall matrix:
            my.matrix <- rbind(my.matrix,my.matrix.p)
        }

        #now for the minus strand:
        tss.vec <- start(x[[i]][strand(x[[i]]) == "-"])
        if (length(tss.vec) > 3) {
# ... no point continuing when there are almost no TSS tags
            my.TSSs <- unique(tss.vec)
            my.matrix.m <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

            this.TSS <- tss.vec[1]
            n.TSSs <- 1
            k <- 0
            for (j in 2:length(tss.vec)) {
                if (tss.vec[j] == this.TSS) {
                    n.TSSs <- n.TSSs + 1
                }
                else {
                    k <- k + 1
                    my.matrix.m[k,] <- c(this.seq, this.TSS, "-", n.TSSs)
                    this.TSS <- tss.vec[j]
                    n.TSSs <- 1
                }
            }
            k <- k + 1
            my.matrix.m[k,] <- c(this.seq, this.TSS, "-", n.TSSs)
# ... adding the minus strand matrix of this.seq to the overall matrix:
            my.matrix <- rbind(my.matrix,my.matrix.m)
        }
    }
    colnames(my.matrix) <- c("seq","TSS","strand","nTSSs")
    my.df <- as.data.frame(my.matrix)
    my.df$seq <- as.character(my.df$seq)
    my.df$TSS <- as.numeric(as.character(my.df$TSS))
    my.df$strand <- as.character(my.df$strand)
    my.df$nTSSs <- as.numeric(as.character(my.df$nTSSs))

    if (writeDF==TRUE) {
        write.table(my.df, outfname, quote=FALSE, col.names=TRUE,
                    row.names=FALSE, sep="\t")
        message("\nThe TSS dataset has been written to file ", outfname,
                "\nin your working directory.")
    }

    return(my.df)
}


################################################################################
#' @title tsrCluster
#' @description an internal function that partitions, then clusters TSS data by
#' sequence to create a data frame containing the coordinates of identified TSRs
#' and other associated metrics, including the count (nTSSs), TSR width and Shape
#' Index (SI). tsrCluster is an internal function that is invoked via detTSR(),
#' which in turn is called by the user-level function determineTSR().
#' 
#' @keywords internal
#'
#' @importFrom gtools mixedsort
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
#' @export

tsrCluster <- function(x, minNbrTSSs=3, minDist=20) {
    tss.df <- x
    uni.seq <- unique(tss.df[,1])
    n.seq <- length(uni.seq)
    tss.df.total <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, tsrWidth=NA, shapeIndex=NA)
    TSS.df <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, tsrWidth=NA, shapeIndex=NA)

    for (l in 1:n.seq) { #by sequence
        seq.name <- mixedsort(uni.seq[l], decreasing=FALSE)
        this.tss <- subset(tss.df, seq==seq.name)
        sTSS <- subset(this.tss, this.tss$nTSSs>=minNbrTSSs)

        #... clustering TSS on the plus strand:
        sTSS.p <- subset(sTSS, strand=="+")
        sTSS.p <- as.matrix(sTSS.p)
        my.len <- nrow(sTSS.p)
        TSS.df.p <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, tsrWidth=NA, shapeIndex=NA)
        if (my.len == 0) {
        }
        else if (my.len == 1) {
            my.tss <- as.numeric(sTSS.p[1,2])
            my.count <- as.numeric(sTSS.p[1,4])
            combined.tss <- rbind(my.tss, my.count)
            tss.out <- tssArrayProperties(combined.tss, seq.name, "+")
            TSS.df.p <- rbind(TSS.df.p, tss.out)
        }
        else {
            my.tss <- as.numeric(sTSS.p[1,2])
            my.count <- as.numeric(sTSS.p[1,4])
            for (i in 1:(my.len-1)) {
                tss.1 <- as.numeric(sTSS.p[i,2])
                tss.1.count <- as.numeric(sTSS.p[i,4])
                tss.2 <- as.numeric(sTSS.p[i+1,2])
                tss.2.count <- as.numeric(sTSS.p[i+1,4])
                tss.dist <- abs(tss.2-tss.1)
                if (tss.dist < minDist) {
                    my.tss <- c(my.tss,tss.2)
                    my.count <- c(my.count, tss.2.count)
                    if (i == my.len-1) {        # wrapping up the last TSR
                        combined.tss <- rbind(my.tss, my.count)
                        tss.out <- tssArrayProperties(combined.tss, seq.name, "+")
                        TSS.df.p <- rbind(TSS.df.p, tss.out)
                    }
                    next
                }
                else {
                    combined.tss <- rbind(my.tss, my.count)
                    tss.out <- tssArrayProperties(combined.tss, seq.name, "+")
                    TSS.df.p <- rbind(TSS.df.p, tss.out)
                    my.tss <- tss.2
                    my.count <- tss.2.count
                }
            }
        }

        #... clustering TSS on the minus strand:

        sTSS.m <- subset(sTSS, strand=="-")
        sTSS.m <- as.matrix(sTSS.m)
        my.len <- nrow(sTSS.m)
        TSS.df.m <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, tsrWidth=NA, shapeIndex=NA)
        if (my.len == 0) {
        }
        else if (my.len == 1) {
            my.tss <- as.numeric(sTSS.m[1,2])
            my.count <- as.numeric(sTSS.m[1,4])
            combined.tss <- rbind(my.tss, my.count)
            tss.out <- tssArrayProperties(combined.tss, seq.name ,"-")
            TSS.df.m <- rbind(TSS.df.m, tss.out)
        }
        else {
            my.tss <- as.numeric(sTSS.m[1,2])
            my.count <- as.numeric(sTSS.m[1,4])
            for (i in 1:(my.len-1)) {
                tss.1 <- as.numeric(sTSS.m[i,2])
                tss.1.count <- as.numeric(sTSS.m[i,4])
                tss.2 <- as.numeric(sTSS.m[i+1,2])
                tss.2.count <- as.numeric(sTSS.m[i+1,4])
                tss.dist <- abs(tss.2-tss.1)
                if (tss.dist < minDist) {
                    my.tss <- c(my.tss,tss.2)
                    my.count <- c(my.count, tss.2.count)
                    if (i == my.len-1) {        # wrapping up the last TSR
                        combined.tss <- rbind(my.tss, my.count)
                        tss.out <- tssArrayProperties(combined.tss, seq.name, "-")
                        TSS.df.m <- rbind(TSS.df.m, tss.out)
                    }
                    next
                }
                else {
                    combined.tss <- rbind(my.tss, my.count)
                    tss.out <- tssArrayProperties(combined.tss, seq.name, "-")
                    TSS.df.m <- rbind(TSS.df.m, tss.out)
                    my.tss <- tss.2
                    my.count <- tss.2.count
                }
            }
        }
        TSS.df.p <- TSS.df.p[-1,]
        TSS.df.m <- TSS.df.m[-1,]
        TSS.df <- rbind(TSS.df.p, TSS.df.m)
        tss.df.total <- rbind(tss.df.total, TSS.df)

    }        #end of by sequence for-loop

    tss.df.total <- tss.df.total[-1,]
    return(tss.df.total)
}


#' @title tssArrayProperties
#' @description An internal function that caculates various properties
#' for a TSR derived in tsrCluster()
#'
#' @keywords internal
#'
#' @return A vector describing the TSR

tssArrayProperties <- function(tssArray, seqName, strand) {
    my.range <- range(tssArray[1,])
    my.counts <- sum(tssArray[2,])
    my.width <- (my.range[2]-my.range[1])+1
    my.SI <- round(TSRshapeIndex(tssArray), digits=2)
    return(c(seqName,my.range[1],my.range[2],strand,my.counts,my.width,my.SI))
}


#' @title TSRshapeIndex
#' @description An internal function that caculates the shape index (SI)
#' of a given TSR from the output of tsrCluster
#'
#' @keywords internal
#'
#' @return Calculates the shape index (SI) for a given TSR

TSRshapeIndex <- function(tssArray) {
    tagcount <- sum(tssArray[2,])
    v <- apply(tssArray, 2,
               function(c) {p <- c[2]/tagcount;
                            return(p * log2(p))
                           }
              )
    SI <- 2 + sum(v)
    return(SI)
}
