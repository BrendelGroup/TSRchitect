#' tssChr (internal function)
#' Retreives tss data from a given experiment by chromosome.
#' @param tssObj an object of class GRanges containing data from a slot of tssData
#' @param chrName the name of the chromosome to select
#' @return a list object containing TSS data for a single chromosome (plus and minus strands)
#' @importFrom BiocGenerics strand start
#' @importFrom GenomeInfoDb seqnames
#' @export

setGeneric(
           name="tssChr",
           def=function(tssObj, chrName) {
               standardGeneric("tssChr")
    }
    )

setMethod("tssChr",
          signature(tssObj="GRanges", chrName="character"),
          function(tssObj, chrName) {
              uni.chr <- as.character(unique(seqnames(tssObj)))
              chr.list <- new("list", plus=numeric(0), minus=numeric(0))
              match_string <- match(chrName, uni.chr)
              if (is.na(match_string)) {
                  stop("The chromosome you selected doesn't exist.")
              }
              this.chr <- uni.chr[uni.chr==chrName]
              #message("\n Extracting tss data from ", this.chr, ".")
              tss.total <- tssObj[seqnames(tssObj)==this.chr,]
              # starting with the TSSs on the plus strand
              tss.plus <- tss.total[as.character(strand(tss.total))=="+",]
              tss.plus.vec <- start(tss.plus)
              chr.list$plus <- tss.plus.vec
              # now for the TSSs on the minus strand
              tss.minus <- tss.total[as.character(strand(tss.total))=="-",]
              tss.minus.vec <- start(tss.minus)
              chr.list$minus <- tss.minus.vec
              tss.total <- NULL
              return(chr.list)
          }
          )

###############################################################################################

#' acquireTSS
#' Retrieves all tss data from a given TSS experiment.
#' @param expName an object of class tssExp containing the tss data
#' @param tssNum the slot number of the tss data to be clustered
#' @return a list object containing TSS data from the entire dataset
#' @importFrom gtools mixedsort
#' @export

setGeneric(
           name="acquireTSS",
           def=function(expName, tssNum) {
               standardGeneric("acquireTSS")
    }
    )

setMethod("acquireTSS",
          signature(expName="tssExp", tssNum="numeric"),
          function(expName, tssNum) {
              cat("\nAcquiring TSS data from sample tssNum = ", tssNum, " ...\n")
              tss.obj <- expName@tssData[[tssNum]]
              uni.chr <- as.character(unique(seqnames(tss.obj)))
              uni.chr <- mixedsort(uni.chr)
              n.chr <- length(uni.chr)
              tss.list <- new("list")
              oldtime <- Sys.time()

              for (i in 1:n.chr) {
                  as.character(uni.chr[i]) -> chrName
                  tssChr(tssObj=tss.obj, chrName) -> tss.out
                  tss.out -> tss.list[[chrName]]
              }
              names(tss.list) <- uni.chr
              newtime <- Sys.time()
              elapsedtime <- newtime - oldtime
              cat("Done with sample tssNum = ", tssNum, " after time: ",print(elapsedtime),".\n\n")
              return(tss.list)
         }
         )

###############################################################################################
#' expressionCTSS
#' Returns a matrix [a, h] where a = the number of unique TSSs and h = the # of tags observed at that position
#' @importFrom gtools mixedsort
#' @export

expressionCTSS <- function(x, dfName="CTSS.txt", writeDF=TRUE) {
        n.chr <- length(names(x)) # how many chromosomes are there in the TSS list?
        uni.chr <- unique(names(x))
        uni.chr <- mixedsort(uni.chr)

        my.matrix <- NULL
        for (i in 1:n.chr) {
#VB Note: Print a progress note on every 10th sequence; 10 should be a parameter
            if (i%%10 == 0) {
                cat("\n... expressionCTSS running with sequence ", i, " of ", n.chr, "\n")
            }
            uni.chr[i] -> this.chr

            #starting with the plus strand:
            tss.vec <- x[[i]]$plus
            if (length(tss.vec) > 3) {	# ... no point continuing when there are almost no TSS tags
                my.CTSSs <- unique(tss.vec)
                my.matrix.p <- matrix(NA, nrow=(length(my.CTSSs)), ncol=4)

                tss.vec[1] -> this.TSS
                1 -> n.TSSs
                0 -> k
                for (j in 2:length(tss.vec)) {
                    if (tss.vec[j] == this.TSS) {
                        n.TSSs + 1 -> n.TSSs
                    }
                    else {
                        k + 1 -> k
                        c(this.chr, this.TSS, n.TSSs,"+") -> my.matrix.p[k,]
                        tss.vec[j] -> this.TSS
                        1 -> n.TSSs
                    }
                }
                k + 1 -> k
                c(this.chr, this.TSS, n.TSSs,"+") -> my.matrix.p[k,]
                my.matrix <- rbind(my.matrix,my.matrix.p)	#adding the plus strand matrix of this.chr to the overall matrix
	    }
            #now for the minus strand:
            if (length(tss.vec) > 3) {	# ... no point continuing when there are almost no TSS tags
                tss.vec <- x[[i]]$minus
                my.CTSSs <- unique(tss.vec)
                my.matrix.m <- matrix(NA, nrow=(length(my.CTSSs)), ncol=4)

                tss.vec[1] -> this.TSS
                1 -> n.TSSs
                0 -> k
                for (j in 2:length(tss.vec)) {
                    if (tss.vec[j] == this.TSS) {
                        n.TSSs + 1 -> n.TSSs
                    }
                    else {
                        k + 1 -> k
                        c(this.chr, this.TSS, n.TSSs,"-") -> my.matrix.m[k,]
                        tss.vec[j] -> this.TSS
                        1 -> n.TSSs
                    }
                }
                k + 1 -> k
                c(this.chr, this.TSS, n.TSSs,"-") -> my.matrix.m[k,]
                my.matrix <- rbind(my.matrix,my.matrix.m)	#adding the minus strand matrix of this.chr to the overall matrix
            }
        }
        colnames(my.matrix) <- c("chr","CTSS","nTSSs","strand")
        my.df <- as.data.frame(my.matrix)
        my.df$CTSS <- as.numeric(as.character(my.df$CTSS))
        my.df$nTSSs <- as.numeric(as.character(my.df$nTSSs))

        if (writeDF==TRUE) {
            write.table(my.df, dfName, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
            message("\nThe TSS dataset has been written to file ", dfName, "\nin your working directory.")
        }

        return(my.df)
        }


##############################################################################################
#' tsrCluster
#' Partitions, then clusters tss data by chromosome. (Internal function)
#' returns a list of TSRs from the data.frame generated by expressionCTSS()
#' @export

.tsrCluster <- function(x, expThresh=5, minDist=20) {
     ctss.df <- x
     ctss.df[,1] <- as.character(ctss.df[,1])
     ctss.df[,4] <- as.character(ctss.df[,4])
     uni.chr <- unique(ctss.df[,1])
     n.chr <- length(uni.chr)
     overall.list <- vector(mode="list", length=n.chr)

     for (l in 1:n.chr) { #by chromosome
         subset(ctss.df, chr==uni.chr[l]) -> this.ctss
         subset(this.ctss, nTSSs>=expThresh) -> sCTSS

         #... clustering TSS on the plus strand:

         subset(sCTSS, strand=="+") -> sCTSS.p
         as.matrix(sCTSS.p) -> sCTSS.p
#VB Note: The following fails if there is a single sCTSS; unclear what the statement
#         is meant to provide - should we ever have NAs in the rows?
#        sCTSS.p[complete.cases(sCTSS.p),] -> sCTSS.p
         nrow(sCTSS.p) -> my.len
	 if (my.len == 0) {
         }
	 else if (my.len == 1) {
             vector(mode="list") -> ctss.list.p
             vector(mode="list") -> counts.list.p
             as.numeric(sCTSS.p[1,2]) -> my.ctss
             as.numeric(sCTSS.p[1,3]) -> my.count
             rbind(my.ctss, my.count) -> combined.ctss
             c("coordinate","count") -> rownames(combined.ctss)
             (1:ncol(combined.ctss)) -> colnames(combined.ctss)
             combined.ctss -> ctss.list.p[[1]]
         }
	 else {
             vector(mode="list") -> ctss.list.p
             vector(mode="list") -> counts.list.p
             as.numeric(sCTSS.p[1,2]) -> my.ctss
             as.numeric(sCTSS.p[1,3]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sCTSS.p[i,2]) -> ctss.1
                 as.numeric(sCTSS.p[i,3]) -> ctss.1.count
                 as.numeric(sCTSS.p[i+1,2]) -> ctss.2
                 as.numeric(sCTSS.p[i+1,3]) -> ctss.2.count
                 abs(ctss.2-ctss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.ctss,ctss.2) -> my.ctss
                     c(my.count, ctss.2.count) -> my.count
                     if (i == my.len-1) {	# wrapping up the last TSR
                         j + 1 -> j
                         rbind(my.ctss, my.count) -> combined.ctss
                         c("coordinate","count") -> rownames(combined.ctss)
                         (1:ncol(combined.ctss)) -> colnames(combined.ctss)
                         combined.ctss -> ctss.list.p[[j]]
                     }
                     next
                 }
                 else {
                     j + 1 -> j
                     rbind(my.ctss, my.count) -> combined.ctss
                     c("coordinate","count") -> rownames(combined.ctss)
                     (1:ncol(combined.ctss)) -> colnames(combined.ctss)
                     combined.ctss -> ctss.list.p[[j]]
                     ctss.2 -> my.ctss
                     ctss.2.count -> my.count
                 }
             }
	 }
         names.len <- length(ctss.list.p)
         names.vec <- vector(mode="character",length=names.len)
         for (k in 1:names.len) {
             paste("tsr", k, sep="") -> names.vec[k]
         }
         names.vec -> names(ctss.list.p)

         #... clustering TSS on the minus strand:

         subset(sCTSS, strand=="-") -> sCTSS.m
         as.matrix(sCTSS.m) -> sCTSS.m
#VB Note: The following fails if there is a single sCTSS; unclear what the statement
#         is meant to provide - should we ever have NAs in the rows?
#        sCTSS.m[complete.cases(sCTSS.m),] -> sCTSS.m
         nrow(sCTSS.m) -> my.len
	 if (my.len == 0) {
         }
	 else if (my.len == 1) {
             vector(mode="list") -> ctss.list.m
             vector(mode="list") -> counts.list.m
             as.numeric(sCTSS.m[1,2]) -> my.ctss
             as.numeric(sCTSS.m[1,3]) -> my.count
             rbind(my.ctss, my.count) -> combined.ctss
             c("coordinate","count") -> rownames(combined.ctss)
             (1:ncol(combined.ctss)) -> colnames(combined.ctss)
             combined.ctss -> ctss.list.m[[1]]
	 }
	 else {
             vector(mode="list") -> ctss.list.m
             vector(mode="list") -> counts.list.m
             as.numeric(sCTSS.m[1,2]) -> my.ctss
             as.numeric(sCTSS.m[1,3]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sCTSS.m[i,2]) -> ctss.1
                 as.numeric(sCTSS.m[i,3]) -> ctss.1.count
                 as.numeric(sCTSS.m[i+1,2]) -> ctss.2
                 as.numeric(sCTSS.m[i+1,3]) -> ctss.2.count
                 abs(ctss.2-ctss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.ctss,ctss.2) -> my.ctss
                     c(my.count, ctss.2.count) -> my.count
                     if (i == my.len-1) {	# wrapping up the last TSR
                         j + 1 -> j
                         rbind(my.ctss, my.count) -> combined.ctss
                         c("coordinate","count") -> rownames(combined.ctss)
                         (1:ncol(combined.ctss)) -> colnames(combined.ctss)
                         combined.ctss -> ctss.list.m[[j]]
                     }
                     next
                 }
                 else {
                     j + 1 -> j
                     rbind(my.ctss, my.count) -> combined.ctss
                     c("coordinate","count") -> rownames(combined.ctss)
                     (1:ncol(combined.ctss)) -> colnames(combined.ctss)
                     combined.ctss -> ctss.list.m[[j]]
                     ctss.2 -> my.ctss
                     ctss.2.count -> my.count
                 }
             }
	 }
         length(ctss.list.m) -> names.len
         vector(mode="character",length=names.len) -> names.vec
         for (k in 1:names.len) {
             paste("tsr", k, sep="") -> names.vec[k]
         }
         names.vec -> names(ctss.list.m)
         list(plus=ctss.list.p, minus=ctss.list.m) -> ctss.list
         ctss.list -> overall.list[[l]]

     }	#end of by chromosome for-loop

     names(overall.list) <- uni.chr
     return(overall.list)
     }
