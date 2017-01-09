#' tssChr (internal function)
#' Retreives tss data from a given experiment by sequence.
#' @param tssObj an object of class GRanges containing data from a slot of tssTagData
#' @param seqName the name of the sequence to select
#' @return a list object containing TSS data for a single sequence (plus and minus strands)
#' @importFrom BiocGenerics strand start
#' @importFrom GenomeInfoDb seqnames
#' @export

setGeneric(
           name="tssChr",
           def=function(tssObj, seqName) {
               standardGeneric("tssChr")
    }
    )

setMethod("tssChr",
          signature(tssObj="GRanges", seqName="character"),
          function(tssObj, seqName) {
              uni.seq <- as.character(unique(seqnames(tssObj)))
              seq.list <- new("list", plus=numeric(0), minus=numeric(0))
              match_string <- match(seqName, uni.seq)
              if (is.na(match_string)) {
                  stop("The sequence you selected doesn't exist.")
              }
              this.seq <- uni.seq[uni.seq==seqName]
              #message("\n Extracting tss data from ", this.seq, ".")
              tss.total <- tssObj[seqnames(tssObj)==this.seq,]
              # starting with the TSSs on the plus strand
              tss.plus <- tss.total[as.character(strand(tss.total))=="+",]
              tss.plus.vec <- start(tss.plus)
              seq.list$plus <- tss.plus.vec
              # now for the TSSs on the minus strand
              tss.minus <- tss.total[as.character(strand(tss.total))=="-",]
              tss.minus.vec <- start(tss.minus)
              seq.list$minus <- tss.minus.vec
              tss.total <- NULL
              return(seq.list)
          }
          )

###############################################################################################

#' acquireTSS
#' Retrieves all tss data from a given TSS experiment.
#' @param experimentName an object of class tssObject containing the tss data
#' @param tssSet the slot number of the tss data to be clustered
#' @return a list object containing TSS data from the entire dataset
#' @importFrom gtools mixedsort
#' @export

setGeneric(
           name="acquireTSS",
           def=function(experimentName, tssSet) {
               standardGeneric("acquireTSS")
    }
    )

setMethod("acquireTSS",
          signature(experimentName="tssObject", tssSet="numeric"),
          function(experimentName, tssSet) {
              cat("\nAcquiring TSS data from sample tssSet = ", tssSet, " ...\n")
              tss.obj <- experimentName@tssTagData[[tssSet]]
              uni.seq <- as.character(unique(seqnames(tss.obj)))
              uni.seq <- mixedsort(uni.seq)
              n.seq <- length(uni.seq)
              tss.list <- new("list")
              oldtime <- Sys.time()

              for (i in 1:n.seq) {
                  as.character(uni.seq[i]) -> seqName
                  tssChr(tssObj=tss.obj, seqName) -> tss.out
                  tss.out -> tss.list[[seqName]]
              }
              names(tss.list) <- uni.seq
              newtime <- Sys.time()
              elapsedtime <- newtime - oldtime
              cat("Done with sample tssSet = ", tssSet, " after time: ",print(elapsedtime),".\n\n")
              return(tss.list)
         }
         )

###############################################################################################
#' tagCountTSS
#' Returns a matrix [a, h] where a = the number of unique TSSs and h = the # of tags observed at that position
#' @importFrom gtools mixedsort
#' @export

tagCountTSS <- function(x, dfName="TSS.txt", writeDF=FALSE) {
        n.seq <- length(names(x)) # how many sequences are there in the TSS list?
        uni.seq <- unique(names(x))
        uni.seq <- mixedsort(uni.seq)

        my.matrix <- NULL
        for (i in 1:n.seq) {
#VB Note: Print a progress note on every 20th sequence; 20 should be a parameter
            if (i%%20 == 0) {
                cat("... tagCountTSS running with sequence ", i, " of ", n.seq, " for TSS set ", dfName, "\n")
            }
            uni.seq[i] -> this.seq

            #starting with the plus strand:
            tss.vec <- x[[i]]$plus
            if (length(tss.vec) > 3) {	# ... no point continuing when there are almost no TSS tags
                my.TSSs <- unique(tss.vec)
                my.matrix.p <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

                tss.vec[1] -> this.TSS
                1 -> n.TSSs
                0 -> k
                for (j in 2:length(tss.vec)) {
                    if (tss.vec[j] == this.TSS) {
                        n.TSSs + 1 -> n.TSSs
                    }
                    else {
                        k + 1 -> k
                        c(this.seq, this.TSS, "+", n.TSSs) -> my.matrix.p[k,]
                        tss.vec[j] -> this.TSS
                        1 -> n.TSSs
                    }
                }
                k + 1 -> k
                c(this.seq, this.TSS, "+", n.TSSs) -> my.matrix.p[k,]
                my.matrix <- rbind(my.matrix,my.matrix.p)	#adding the plus strand matrix of this.seq to the overall matrix
            }

            #now for the minus strand:
            tss.vec <- x[[i]]$minus
            if (length(tss.vec) > 3) {	# ... no point continuing when there are almost no TSS tags
                my.TSSs <- unique(tss.vec)
                my.matrix.m <- matrix(NA, nrow=(length(my.TSSs)), ncol=4)

                tss.vec[1] -> this.TSS
                1 -> n.TSSs
                0 -> k
                for (j in 2:length(tss.vec)) {
                    if (tss.vec[j] == this.TSS) {
                        n.TSSs + 1 -> n.TSSs
                    }
                    else {
                        k + 1 -> k
                        c(this.seq, this.TSS, "-", n.TSSs) -> my.matrix.m[k,]
                        tss.vec[j] -> this.TSS
                        1 -> n.TSSs
                    }
                }
                k + 1 -> k
                c(this.seq, this.TSS, "-", n.TSSs) -> my.matrix.m[k,]
                my.matrix <- rbind(my.matrix,my.matrix.m)	#adding the minus strand matrix of this.seq to the overall matrix
            }
        }
        colnames(my.matrix) <- c("seq","TSS","strand","nTSSs")
        my.df <- as.data.frame(my.matrix)
        my.df$seq <- as.character(my.df$seq)
        my.df$TSS <- as.numeric(as.character(my.df$TSS))
        my.df$strand <- as.character(my.df$strand)
        my.df$nTSSs <- as.numeric(as.character(my.df$nTSSs))

        if (writeDF==TRUE) {
            write.table(my.df, dfName, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
            message("\nThe TSS dataset has been written to file ", dfName, "\nin your working directory.")
        }

        return(my.df)
        }


##############################################################################################
#' tsrCluster
#' Partitions, then clusters tss data by sequence. (Internal function)
#' returns a list of TSRs from the data.frame generated by tagCountTSS()
#' @export

tsrCluster <- function(x, minNbrTSSs=3, minDist=20) {
     tss.df <- x
     uni.seq <- unique(tss.df[,1])
     n.seq <- length(uni.seq)
     overall.list <- vector(mode="list", length=n.seq)

     for (l in 1:n.seq) { #by sequence
         subset(tss.df, seq==uni.seq[l]) -> this.tss
         subset(this.tss, nTSSs>=minNbrTSSs) -> sTSS

         #... clustering TSS on the plus strand:

         subset(sTSS, strand=="+") -> sTSS.p
         as.matrix(sTSS.p) -> sTSS.p
         nrow(sTSS.p) -> my.len
         if (my.len == 0) { # ... assign an empty list to tss.list.p
             vector(mode="list") -> tss.list.p
         }
         else if (my.len == 1) {
             vector(mode="list") -> tss.list.p
             as.numeric(sTSS.p[1,2]) -> my.tss
             as.numeric(sTSS.p[1,4]) -> my.count
             rbind(my.tss, my.count) -> combined.tss
             c("coordinate","count") -> rownames(combined.tss)
             (1:ncol(combined.tss)) -> colnames(combined.tss)
             combined.tss -> tss.list.p[[1]]
         }
         else {
             vector(mode="list") -> tss.list.p
             as.numeric(sTSS.p[1,2]) -> my.tss
             as.numeric(sTSS.p[1,4]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sTSS.p[i,2]) -> tss.1
                 as.numeric(sTSS.p[i,4]) -> tss.1.count
                 as.numeric(sTSS.p[i+1,2]) -> tss.2
                 as.numeric(sTSS.p[i+1,4]) -> tss.2.count
                 abs(tss.2-tss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.tss,tss.2) -> my.tss
                     c(my.count, tss.2.count) -> my.count
                     if (i == my.len-1) {	# wrapping up the last TSR
                         j + 1 -> j
                         rbind(my.tss, my.count) -> combined.tss
                         c("coordinate","count") -> rownames(combined.tss)
                         (1:ncol(combined.tss)) -> colnames(combined.tss)
                         combined.tss -> tss.list.p[[j]]
                     }
                     next
                 }
                 else {
                     j + 1 -> j
                     rbind(my.tss, my.count) -> combined.tss
                     c("coordinate","count") -> rownames(combined.tss)
                     (1:ncol(combined.tss)) -> colnames(combined.tss)
                     combined.tss -> tss.list.p[[j]]
                     tss.2 -> my.tss
                     tss.2.count -> my.count
                 }
             }
         }
         if (my.len > 0) {
             names.len <- length(tss.list.p)
             names.vec <- vector(mode="character",length=names.len)
             for (k in 1:names.len) {
                 paste("tsr", k, sep="") -> names.vec[k]
             }
             names.vec -> names(tss.list.p)
         }

         #... clustering TSS on the minus strand:

         subset(sTSS, strand=="-") -> sTSS.m
         as.matrix(sTSS.m) -> sTSS.m
         nrow(sTSS.m) -> my.len
         if (my.len == 0) { # ... assign an empty list to tss.list.m
             vector(mode="list") -> tss.list.m
         }
         else if (my.len == 1) {
             vector(mode="list") -> tss.list.m
             as.numeric(sTSS.m[1,2]) -> my.tss
             as.numeric(sTSS.m[1,4]) -> my.count
             rbind(my.tss, my.count) -> combined.tss
             c("coordinate","count") -> rownames(combined.tss)
             (1:ncol(combined.tss)) -> colnames(combined.tss)
             combined.tss -> tss.list.m[[1]]
         }
         else {
             vector(mode="list") -> tss.list.m
             as.numeric(sTSS.m[1,2]) -> my.tss
             as.numeric(sTSS.m[1,4]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sTSS.m[i,2]) -> tss.1
                 as.numeric(sTSS.m[i,4]) -> tss.1.count
                 as.numeric(sTSS.m[i+1,2]) -> tss.2
                 as.numeric(sTSS.m[i+1,4]) -> tss.2.count
                 abs(tss.2-tss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.tss,tss.2) -> my.tss
                     c(my.count, tss.2.count) -> my.count
                     if (i == my.len-1) {	# wrapping up the last TSR
                         j + 1 -> j
                         rbind(my.tss, my.count) -> combined.tss
                         c("coordinate","count") -> rownames(combined.tss)
                         (1:ncol(combined.tss)) -> colnames(combined.tss)
                         combined.tss -> tss.list.m[[j]]
                     }
                     next
                 }
                 else {
                     j + 1 -> j
                     rbind(my.tss, my.count) -> combined.tss
                     c("coordinate","count") -> rownames(combined.tss)
                     (1:ncol(combined.tss)) -> colnames(combined.tss)
                     combined.tss -> tss.list.m[[j]]
                     tss.2 -> my.tss
                     tss.2.count -> my.count
                 }
             }
         }
         if (my.len > 0) {
             length(tss.list.m) -> names.len
             vector(mode="character",length=names.len) -> names.vec
             for (k in 1:names.len) {
                 paste("tsr", k, sep="") -> names.vec[k]
             }
             names.vec -> names(tss.list.m)
         }
         list(plus=tss.list.p, minus=tss.list.m) -> tss.list
         tss.list -> overall.list[[l]]

     }	#end of by sequence for-loop

     names(overall.list) <- uni.seq
     return(overall.list)
     }
