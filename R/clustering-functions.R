
#' tssChr (internal function)
#' Retreives tss data from a given experiment by chromosome.
#' @param expName an object of class tssExp containing the tss data
#' @param tssNum the slot number of the tss data to be clustered
#' @param chrName the name of the chromosome to select
#' @return a list object containing TSS data for a single chromosome (plus and minus strands)
#' @importFrom BiocGenerics strand start
#' @importFrom GenomeInfoDb seqnames
#' @export

setGeneric(
           name="tssChr",
           def=function(expName, tssNum, chrName) {
               standardGeneric("tssChr")
    }
    )

setMethod("tssChr",
          signature(expName="tssExp", tssNum="numeric", chrName="character"),

          function(expName, tssNum, chrName) {

              if (tssNum>length(expName@tssData)) {

                  stop("The value selected exceeds the number of slots in tssData.")

              }

              this.tss <- expName@tssData[[tssNum]]

              uni.chr <- as.character(unique(seqnames(this.tss)))

              chr.list <- new("list", plus=numeric(0), minus=numeric(0))

              match_string <- match(chrName, uni.chr)
              
              if (is.na(match_string)) {

                  stop("The chromosome you selected doesn't exist.")

              }
              
                  this.chr <- uni.chr[uni.chr==chrName]

#                  message("\n Extracting tss data from ", this.chr, ".")

                  tss.total <- this.tss[seqnames(this.tss)==this.chr,]

                  # starting with the TSSs on the plus strand
                  
                  tss.plus <- tss.total[strand(tss.total)=="+",]

                  tss.plus.vec <- start(tss.plus)

                  chr.list$plus <- tss.plus.vec

                  # now for the TSSs on the minus strand

                  tss.minus <- tss.total[strand(tss.total)=="-",]

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

              this.tss <- expName@tssData[[tssNum]]

              uni.chr <- as.character(unique(seqnames(this.tss)))

              uni.chr <- mixedsort(uni.chr)

              n.chr <- length(uni.chr)

              tss.list <- new("list")

              for (i in 1:n.chr) {

                  uni.chr[i] -> chrName

                  tssChr(expName, tssNum, chrName) -> tss.out

                  tss.out -> tss.list[[i]]

                  }
                 
            names(tss.list) <- uni.chr

            return(tss.list)
          }
           )   

###############################################################################################
#' expressionCTSS
#' Returns a matrix [a, h] where a = the number of unique TSSs and h = the # of tags observed at that position
#' @importFrom gtools mixedsort
#' @export

expressionCTSS <- function(x, writeDF=TRUE, dfName="CTSS.txt") {

        #starting with the plus strand

        my.matrix <- matrix(NA, nrow=1, ncol=4)
        n.chr <- length(names(x)) # how many chromosomes are there in the TSS list?
        uni.chr <- unique(names(x))
        uni.chr <- mixedsort(uni.chr)
        
        for (i in 1:n.chr) {

            uni.chr[i] -> this.chr

#            ptm <- proc.time() #for timing
#            cat("Transforming the data from TSS dataset", this.chr,"\n")

            tss.vec <- x[[i]]$plus
            my.CTSSs <- unique(tss.vec)
            my.matrix.p <- matrix(NA, nrow=(length(my.CTSSs)), ncol=4)
            for (j in 1:length(my.CTSSs)) {
                my.CTSSs[j] -> this.TSS
                which(tss.vec==this.TSS) -> my.ind
                length(my.ind) -> n.TSSs
                c(this.chr, this.TSS, n.TSSs,"+") -> my.matrix.p[j,]
            }

            #now for the minus strand

            tss.vec <- x[[i]]$minus
            my.CTSSs <- unique(tss.vec)
            my.CTSSs <- sort(my.CTSSs)
            my.matrix.m <- matrix(NA, nrow=(length(my.CTSSs)), ncol=4)
            
            for (j in 1:length(my.CTSSs)) {
                my.CTSSs[j] -> this.TSS
                which(tss.vec==this.TSS) -> my.ind
                length(my.ind) -> n.TSSs
                c(this.chr, this.TSS, n.TSSs, "-") -> my.matrix.m[j,]
            }
#            print(proc.time() - ptm) #for timing

            this.matrix <- rbind(my.matrix.p, my.matrix.m) #combining the two matrices
            my.matrix <- rbind(my.matrix, this.matrix)
        }

        colnames(my.matrix) <- c("chr","CTSS","nTSSs","strand")
        my.matrix <- my.matrix[-1,] #removing the first row, which contains only NAs
        my.df <- as.data.frame(my.matrix)
        my.df$CTSS <- as.numeric(as.character(my.df$CTSS))
        my.df$nTSSs <- as.numeric(as.character(my.df$nTSSs))

        if (writeDF==TRUE) {
            write.table(my.df, dfName, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
        }
        return(my.df)
        }
            
##############################################################################################
#' tsrCluster
#' Partitions, then clusters tss data by chromosome. (Internal function)
#' @export

.tsrCluster <- function(x, expThresh=5, minDist=20) {
    ## returns a list of TSRs from the data.frame generated by expressionCTSS()
    
     ctss.df <- x
     ctss.df[,1] <- as.character(ctss.df[,1])
     ctss.df[,2] <- as.numeric(as.character(ctss.df[,2]))
     ctss.df[,3] <- as.numeric(as.character(ctss.df[,3]))
     ctss.df[,4] <- as.character(ctss.df[,4])
     chr.vec <- as.character(ctss.df[,1])
     uni.chr <- unique(chr.vec)
     n.chr <- length(uni.chr)
     overall.list <- vector(mode="list")

     for (l in 1:n.chr) { #by chromosome
         as.character(uni.chr[l]) -> my.chr
         subset(ctss.df, chr== my.chr) -> this.ctss #starting with plus
         subset(ctss.df, nTSSs>= expThresh) -> sCTSS
         subset(sCTSS, strand>= "+") -> sCTSS.p
         as.matrix(sCTSS.p) -> sCTSS.p #a kludge we'll use for now
         sCTSS.p[complete.cases(sCTSS.p),] -> sCTSS.p
         nrow(sCTSS.p) -> my.len
         vector(mode="list") -> ctss.list.p
         as.numeric(sCTSS.p[1,2]) -> my.ctss
         0 -> j
         for (i in 1:(my.len-1)) {
            as.numeric(sCTSS.p[i,2]) -> ctss.1 
            as.numeric(sCTSS.p[(i+1),2]) -> ctss.2
            abs(ctss.2-ctss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.ctss,ctss.2) -> my.ctss
                     next
                 }
                 else {
                     j + 1 -> j 
                     my.ctss -> ctss.list.p[[j]]
                     ctss.2 -> my.ctss
                }
        }
     names.len <- length(ctss.list.p)
     names.vec <- vector(mode="character",length=names.len)
     for (k in 1:names.len) {
         paste("tsr", k, sep="") -> names.vec[k]
     }
     names(ctss.list.p) <- names.vec

     sCTSS.m <- subset(sCTSS, strand>= "-") #now with the minus strand
     sCTSS.m <- as.matrix(sCTSS.m) #a kludge we'll use for now
     sCTSS.m <- sCTSS.m[complete.cases(sCTSS.m),] 
     my.len <- nrow(sCTSS.m)
     ctss.list.m <- vector(mode="list")
     as.numeric(sCTSS.m[1,2]) -> my.ctss
     0 -> j
     for (i in 1:(my.len-1)) {
            as.numeric(sCTSS.m[i,2]) -> ctss.1 
            as.numeric(sCTSS.m[(i+1),2]) -> ctss.2
            abs(ctss.2-ctss.1) -> tss.dist
                 if (tss.dist < minDist) {
                     c(my.ctss,ctss.2) -> my.ctss
                     next
                 }
                 else {
                     j + 1 -> j 
                     my.ctss -> ctss.list.m[[j]]
                     ctss.2 -> my.ctss
                }
        }
     names.len <- length(ctss.list.m)
     names.vec <- vector(mode="character",length=names.len)
     for (k in 1:names.len) {
         paste("tsr", k, sep="") -> names.vec[k]
     }
     names(ctss.list.m) <- names.vec
     ctss.list <- list(plus=ctss.list.p, minus=ctss.list.m)
     overall.list[[l]] <- ctss.list
         
     }
     names(overall.list) <- uni.chr
     return(overall.list)
 }

