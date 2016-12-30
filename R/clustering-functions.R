#' tssChr (internal function)
#' Retreives tss data from a given experiment by chromosome.
#' @param tssObj an object of class GRanges containing data from a slot of tssTagData
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
        n.chr <- length(names(x)) # how many chromosomes are there in the TSS list?
        uni.chr <- unique(names(x))
        uni.chr <- mixedsort(uni.chr)

        my.matrix <- NULL
        for (i in 1:n.chr) {
#VB Note: Print a progress note on every 20th sequence; 20 should be a parameter
            if (i%%20 == 0) {
                cat("... tagCountTSS running with sequence ", i, " of ", n.chr, " for TSS set ", dfName, "\n")
            }
            uni.chr[i] -> this.chr

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
        colnames(my.matrix) <- c("chr","TSS","nTSSs","strand")
        my.df <- as.data.frame(my.matrix)
        my.df$chr <- as.character(my.df$chr)
        my.df$TSS <- as.numeric(as.character(my.df$TSS))
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
#' returns a list of TSRs from the data.frame generated by tagCountTSS()
#' @export

tsrCluster <- function(x, minNbrTSSs=3, minDist=20) {
     tss.df <- x
     tss.df[,1] <- as.character(tss.df[,1])
     tss.df[,4] <- as.character(tss.df[,4])
     uni.chr <- unique(tss.df[,1])
     n.chr <- length(uni.chr)
     overall.list <- vector(mode="list", length=n.chr)

     for (l in 1:n.chr) { #by chromosome
         subset(tss.df, chr==uni.chr[l]) -> this.tss
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
             as.numeric(sTSS.p[1,3]) -> my.count
             rbind(my.tss, my.count) -> combined.tss
             c("coordinate","count") -> rownames(combined.tss)
             (1:ncol(combined.tss)) -> colnames(combined.tss)
             combined.tss -> tss.list.p[[1]]
         }
	 else {
             vector(mode="list") -> tss.list.p
             as.numeric(sTSS.p[1,2]) -> my.tss
             as.numeric(sTSS.p[1,3]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sTSS.p[i,2]) -> tss.1
                 as.numeric(sTSS.p[i,3]) -> tss.1.count
                 as.numeric(sTSS.p[i+1,2]) -> tss.2
                 as.numeric(sTSS.p[i+1,3]) -> tss.2.count
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
             as.numeric(sTSS.m[1,3]) -> my.count
             rbind(my.tss, my.count) -> combined.tss
             c("coordinate","count") -> rownames(combined.tss)
             (1:ncol(combined.tss)) -> colnames(combined.tss)
             combined.tss -> tss.list.m[[1]]
	 }
	 else {
             vector(mode="list") -> tss.list.m
             as.numeric(sTSS.m[1,2]) -> my.tss
             as.numeric(sTSS.m[1,3]) -> my.count
             0 -> j
             for (i in 1:(my.len-1)) {
                 as.numeric(sTSS.m[i,2]) -> tss.1
                 as.numeric(sTSS.m[i,3]) -> tss.1.count
                 as.numeric(sTSS.m[i+1,2]) -> tss.2
                 as.numeric(sTSS.m[i+1,3]) -> tss.2.count
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

     }	#end of by chromosome for-loop

     names(overall.list) <- uni.chr
     return(overall.list)
     }
