################################################################################
#' tagCountTSS
#' @description an internal function that returns a matrix [a, h] where a = the
#' number of unique TSSs and h = the # of tags observed at that position.
#' tagCountTSS is invoked via prcTSS.
#'
#' @import BiocGenerics
#' @import GenomicRanges
#' @importFrom utils write.table
#'
#' @param y a data frame containing the contents of a single slot of tssTagData.
#' @param n.cores the number of cores to be used for this job.
#' @param outfname the prefix of the file name of TSS information to be written.
#' (character)
#' @param writeTable if TRUE, the tag count information is written to a file in the
#' workding directory (logical)
#'
#' @keywords internal
#' @return a matrix [a, h] containing the number of unique TSSs (a) and
#' their corresponding abundances (h) is returned.


tagCountTSS <- function(y, n.cores=1, outfname="TSS.txt", writeTable=FALSE) {
    x <- S4Vectors::split(y,seqnames(y))
    combined.matrix <- NULL

    fctz <- function(z) {
        my.matrix <- NULL
        this.seq <- as.character(z@seqnames@values[1])
        if (is.na(this.seq)) {		#z is empty (no tags for the sequence)
            return(my.matrix)
	}

        #starting with the plus strand:
        tss.vec <- start(z[strand(z) == "+"])
        if (length(tss.vec) > 3) {        #stop if there are nearly no tags
            my.TSSs <- unique(tss.vec)
            my.matrix.p <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

            this.TSS <- tss.vec[1]
            n.TAGs <- 1
            k <- 0
            for (j in 2:length(tss.vec)) {
                if (tss.vec[j] == this.TSS) {
                    n.TAGs <- n.TAGs + 1
                }
                else {
                    k <- k + 1
                    my.matrix.p[k,] <- c(this.seq, this.TSS, "+", n.TAGs)
                    this.TSS <- tss.vec[j]
                    n.TAGs <- 1
                }
            }
            k <- k + 1
            my.matrix.p[k,] <- c(this.seq, this.TSS, "+", n.TAGs)
# ... add the plus strand matrix of this.seq to the overall matrix:
            my.matrix <- rbind(my.matrix,my.matrix.p)
        }

        #now for the minus strand:
        tss.vec <- start(z[strand(z) == "-"])
        if (length(tss.vec) > 3) {
# ... no point continuing when there are almost no TSS tags
            my.TSSs <- unique(tss.vec)
            my.matrix.m <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

            this.TSS <- tss.vec[1]
            n.TAGs <- 1
            k <- 0
            for (j in 2:length(tss.vec)) {
                if (tss.vec[j] == this.TSS) {
                    n.TAGs <- n.TAGs + 1
                }
                else {
                    k <- k + 1
                    my.matrix.m[k,] <- c(this.seq, this.TSS, "-", n.TAGs)
                    this.TSS <- tss.vec[j]
                    n.TAGs <- 1
                }
            }
            k <- k + 1
            my.matrix.m[k,] <- c(this.seq, this.TSS, "-", n.TAGs)
# ... adding the minus strand matrix of this.seq to the overall matrix:
            my.matrix <- rbind(my.matrix,my.matrix.m)
        }
	return(my.matrix)
    }

    lm <- lapply(x,fctz)
    combined.matrix <- do.call(rbind,lm)
    colnames(combined.matrix) <- c("seq","TSS","strand","nTAGs")
    my.df <- as.data.frame(combined.matrix)
    my.df$seq <- as.character(my.df$seq)
    my.df$TSS <- as.numeric(as.character(my.df$TSS))
    my.df$strand <- as.character(my.df$strand)
    my.df$nTAGs <- as.numeric(as.character(my.df$nTAGs))

    if (writeTable==TRUE) {
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
#' and other associated metrics, including the site and tag counts (nTSSs, nTAGs),
#' TSR width and (modified) Shape Index (SI, mSI). tsrCluster is an internal
#' function that is invoked via detTSR(), which in turn is called by the
#' user-level function determineTSR().
#' 
#' @keywords internal
#'
#' @importFrom gtools mixedsort
#'
#' @param x a data frame containing a single slot from either tssCountData or
#' tssCountDataMerged, depending on its invocation by parent function
#' determineTSR()
#' @param minNbrTAGs the minimum number of tags at a given TSS position
#' for a TSS to be included in clustering. (numeric)
#' @param minDist the maximum distance of TSSs between two TSRs in base pairs.
#' (numeric)
#'
#' @return A data frame of TSRs with variables\cr
#' \enumerate{
#'         \item seq       = sequence identifier (seq)
#'         \item start     = start of TSR (num)
#'         \item end       = end of TSR (num)
#'         \item strand    = + or - (factor)
#'         \item nTSSs     = count of TSSs (num)
#'         \item nTAGs     = count of tags (num)
#'         \item tsrPeak   = maximum tag count fraction for all TSS positions
#'                           in the TSR (num)
#'         \item tsrWdth   = width of TSR (num)
#'         \item tsrTrq    = TSR torque; measure of TSR balance (num)
#'         \item tsrSI     = shape index value of TSR (num)
#'         \item tsrMSI    = modified shape index value of TSR (num)
#' }
#' @export

tsrCluster <- function(x, minNbrTAGs=3, minDist=20) {
    tss.df <- x
    uni.seq <- mixedsort(unique(tss.df[,1]), decreasing=FALSE)
    tss.df.total <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, nTAGs=NA,
			       tsrPeak=NA, tsrWdth=NA, tsrTrq=NA,
			       tsrSI=NA, tsrMSI=NA)
    TSR.df <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, nTAGs=NA,
			       tsrPeak=NA, tsrWdth=NA, tsrTrq=NA,
			       tsrSI=NA, tsrMSI=NA)

    fctn <- function(seqname) {
        this.tss <- subset(tss.df, seq==seqname)
        sTSS <- subset(this.tss, this.tss$nTAGs>=minNbrTAGs)

        #... clustering TSS on the plus strand:
        sTSS.p <- subset(sTSS, strand=="+")
        sTSS.p <- as.matrix(sTSS.p)
        my.len <- nrow(sTSS.p)
        TSR.df.p <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, nTAGs=NA,
			       tsrPeak=NA, tsrWdth=NA, tsrTrq=NA,
			       tsrSI=NA, tsrMSI=NA)
        if (my.len == 0) {	# ... nothing to do
        }
        else if (my.len == 1) {
            my.tss <- as.numeric(sTSS.p[1,2])
            my.count <- as.numeric(sTSS.p[1,4])
            clustered.tss <- cbind(my.tss, my.count)
            my.tsr <- tssArrayProperties(clustered.tss, seqname, "+")
            TSR.df.p <- rbind(TSR.df.p, my.tsr)
        }
        else {
            my.tss <- as.numeric(sTSS.p[1,2])
            my.count <- as.numeric(sTSS.p[1,4])
            my.nbrtss <- 1
            for (i in 1:(my.len-1)) {
                tss.1 <- as.numeric(sTSS.p[i,2])
                tss.1.count <- as.numeric(sTSS.p[i,4])
                tss.2 <- as.numeric(sTSS.p[i+1,2])
                tss.2.count <- as.numeric(sTSS.p[i+1,4])
                tss.dist <- abs(tss.2-tss.1)
                if (tss.dist < minDist) {
                    if (tss.dist == 0) {
                        my.count[my.nbrtss] <- my.count[my.nbrtss] + tss.2.count
                    }
                    else {
                        my.tss <- rbind(my.tss,tss.2)
                        my.count <- rbind(my.count, tss.2.count)
                        my.nbrtss <- my.nbrtss + 1
                    }
                    if (i == my.len-1) {        # wrapping up the last TSR
                        clustered.tss <- cbind(my.tss, my.count)
                        my.tsr <- tssArrayProperties(clustered.tss, seqname, "+")
                        TSR.df.p <- rbind(TSR.df.p, my.tsr)
                    }
                    next
                }
                else {
                    clustered.tss <- cbind(my.tss, my.count)
                    my.tsr <- tssArrayProperties(clustered.tss, seqname, "+")
                    TSR.df.p <- rbind(TSR.df.p, my.tsr)
                    my.tss <- tss.2
                    my.count <- tss.2.count
                    my.nbrtss <- 1
                }
            }
        }

        #... clustering TSS on the minus strand:

        sTSS.m <- subset(sTSS, strand=="-")
        sTSS.m <- as.matrix(sTSS.m)
        my.len <- nrow(sTSS.m)
        TSR.df.m <- data.frame(seq=NA, start=NA, end=NA, strand=NA,
                               nTSSs=NA, nTAGs=NA, tsrPeak=NA, tsrWdth=NA,
                               tsrTrq=NA, tsrSI=NA, tsrMSI=NA)
        if (my.len == 0) {
        }
        else if (my.len == 1) {
            my.tss <- as.numeric(sTSS.m[1,2])
            my.count <- as.numeric(sTSS.m[1,4])
            clustered.tss <- cbind(my.tss, my.count)
            my.tsr <- tssArrayProperties(clustered.tss, seqname,"-")
            TSR.df.m <- rbind(TSR.df.m, my.tsr)
        }
        else {
            my.tss <- as.numeric(sTSS.m[1,2])
            my.count <- as.numeric(sTSS.m[1,4])
            my.nbrtss <- 1
            for (i in 1:(my.len-1)) {
                tss.1 <- as.numeric(sTSS.m[i,2])
                tss.1.count <- as.numeric(sTSS.m[i,4])
                tss.2 <- as.numeric(sTSS.m[i+1,2])
                tss.2.count <- as.numeric(sTSS.m[i+1,4])
                tss.dist <- abs(tss.2-tss.1)
                if (tss.dist < minDist) {
                    if (tss.dist == 0) {
                        my.count[my.nbrtss] <- my.count[my.nbrtss] + tss.2.count
                    }
                    else {
                        my.tss <- rbind(my.tss,tss.2)
                        my.count <- rbind(my.count, tss.2.count)
                        my.nbrtss <- my.nbrtss + 1
                    }
                    if (i == my.len-1) {        # wrapping up the last TSR
                        clustered.tss <- cbind(my.tss, my.count)
                        my.tsr <- tssArrayProperties(clustered.tss, seqname, "-")
                        TSR.df.m <- rbind(TSR.df.m, my.tsr)
                    }
                    next
                }
                else {
                    clustered.tss <- cbind(my.tss, my.count)
                    my.tsr <- tssArrayProperties(clustered.tss, seqname, "-")
                    TSR.df.m <- rbind(TSR.df.m, my.tsr)
                    my.tss <- tss.2
                    my.count <- tss.2.count
                    my.nbrtss <- 1
                }
            }
        }
        TSR.df.p <- TSR.df.p[-1,]
        TSR.df.m <- TSR.df.m[-1,]
        TSR.df <- rbind(TSR.df.p, TSR.df.m)
        return(TSR.df)
    }

    ldf <- lapply(uni.seq,fctn)
    tss.df.total <- do.call(rbind,ldf)
    return(tss.df.total)
}


#' @title tssArrayProperties
#' @description An internal function that calculates various properties
#' for a TSR derived in tsrCluster()
#'
#' @keywords internal
#'
#' @param tssArray an object containing TSS coordinates and their
#' abundances. (data.frame)
#' @param seqName the name of the chromosome or scaffold. (character)
#' @param strand the strand that the TSR tags are located. (character)
#'
#' @return A vector containing information about the TSR.
#' The returned vector is as follows:
#' seqName (character), TSR start (numeric), TSR end (numeric), strand (character),
#' number of TSSs (numeric), number of tags (numeric), fraction of tags in highest
#' peak (numeric), TSR width (numeric), TSR torque (numeric),
#' Shape Index (numeric), Modified Shape Index (numeric)

tssArrayProperties <- function(tssArray, seqName, strand) {
    tsr.range    <- range(tssArray[,1])
    tsr.midpoint <- (tsr.range[1] + tsr.range[2])/2
    tsr.tsscount <- length(tssArray[,1])
    tsr.tagcount <- sum(tssArray[,2])
    tsr.peak     <- round(max(tssArray[,2])/tsr.tagcount, digits=2)
    tsr.torque   <- round(sum((tssArray[,2]/tsr.tagcount) *
			      (tssArray[,1]-tsr.midpoint) ), digits=2)
    tsr.width    <- (tsr.range[2]-tsr.range[1])+1
    tsr.SI       <- round(TSRshapeIndex(tssArray), digits=2)
    tsr.mSI      <- round(TSRmshapeIndex(tssArray), digits=2)
    return(list(seqName,tsr.range[1],tsr.range[2],strand,tsr.tsscount,
                tsr.tagcount,tsr.peak,tsr.width,tsr.torque,tsr.SI,tsr.mSI))
}


#' @title TSRshapeIndex
#' @description An internal function that caculates the shape index (SI)
#' of a given TSR from the output of tsrCluster.
#' TSRshapeIndex is called by tssArrayProperties().
#'
#' @keywords internal
#'
#' @param tssArray an object containing TSS coordinates and their
#' abundances. (data.frame)
#'
#' @return Returns a  numeric vector of length 1 containing the Shape Index (SI)
#' value for the selected TSR.

TSRshapeIndex <- function(tssArray) {
    tagcount <- sum(tssArray[,2])
    v <- apply(tssArray, 1,
               function(c) {p <- c[2]/tagcount;
                            return(p * log2(p))
                           }
              )
    SI <- 2 + sum(v)
    return(SI)
}


#' @title TSRmshapeIndex
#' @description An internal function that caculates the modified shape index (mSI)
#' of a given TSR from the output of tsrCluster.
#' TSRmshapeIndex is called by tssArrayProperties().
#'
#' @keywords internal
#'
#' @param tssArray an object containing TSS coordinates and their
#' abundances. (data.frame)
#'
#' @return Returns a  numeric vector of length 1 containing the Modified Shape Index (mSI)
#' value for the selected TSR.

TSRmshapeIndex <- function(tssArray) {
    tagcount <- sum(tssArray[,2])
    tsscount <- length(tssArray[,1])
    v <- apply(tssArray, 1,
               function(c) {p <- c[2]/tagcount;
                    	    if (tsscount == 1) {
                              return(0)
	                    } else {
                              return(p * log(p)/log(tsscount))
	                    }
                           }
              )
    mSI <- 1 + sum(v)
    return(mSI)
}
