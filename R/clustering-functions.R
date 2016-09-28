
#' Partitions, then clusters tss data by chromosome.
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

              chr.list <- list(plus="plus", minus="minus")

              match_string <- match(chrName, uni.chr)
              
              if (is.na(match_string)) {

                  stop("The chromosome you selected doesn't exist.")

              }
              
                  this.chr <- uni.chr[uni.chr==chrName]

                  message("\n Extracting tss data from ", this.chr, ".\n")

                  tss.total <- this.tss[seqnames(this.tss)==this.chr,]

                  # starting with the TSSs on the plus strand
                  
                  tss.plus <- tss.total[strand(tss.total)=="+",]

                  tss.plus.vec <- start(tss.plus)

                  chr.list$plus$tss <- tss.plus.vec

                  # now for the TSSs on the minus strand

                  tss.minus <- tss.total[strand(tss.total)=="-",]

                  tss.minus.vec <- start(tss.minus)

                  chr.list$minus$tss <- tss.minus.vec

                  tss.total <- NULL

                  return(chr.list)
          }
           )   

###############################################################################################

expressionCTSS <- function(x) {
        ## returns a matrix [a, h] where a = the number of unique TSSs and h = the # of tags observed at that position

        #starting with the plus strand

#            ptm <- proc.time() #for timing

            tss.vec <- x$plus$tss
            my.CTSSs <- unique(tss.vec)
            my.matrix.p <- matrix(NA, nrow=(length(my.CTSSs)), ncol=3)
            for (i in 1:length(my.CTSSs)) {
                my.CTSSs[i] -> this.TSS
                which(tss.vec==this.TSS) -> my.ind
                length(my.ind) -> n.TSSs
                c(this.TSS, n.TSSs,"+") -> my.matrix.p[i,]
            }

#            print(proc.time() - ptm)
            
            #now for the minus strand

            ptm <- proc.time() #for timing

            tss.vec <- x$minus$tss
            my.CTSSs <- unique(tss.vec)
            my.CTSSs <- sort(my.CTSSs)
            my.matrix.m <- matrix(NA, nrow=(length(my.CTSSs)), ncol=3)
            
            for (i in 1:length(my.CTSSs)) {
                my.CTSSs[i] -> this.TSS
                which(tss.vec==this.TSS) -> my.ind
                length(my.ind) -> n.TSSs
                c(this.TSS, n.TSSs, "-") -> my.matrix.m[i,]
            }

#            print(proc.time() - ptm) #for timing
            
            my.matrix <- rbind(my.matrix.p, my.matrix.m) #combining the two matrices
            colnames(my.matrix) <- c("CTSS","nTSSs","strand")
            my.df <- as.data.frame(my.matrix)

            return(my.df)
        }
            
###############################################################################################

tsrCluster <- function(x, expThresh=5, minDist=20) {
    ## returns a list of TSRs from a given chromosome or scaffold

     ctss.df <- x
     ctss.df[,1] <- as.numeric(as.character(ctss.df[,1]))
     ctss.df[,2] <- as.numeric(as.character(ctss.df[,2]))
     ctss.df[,3] <- as.character(ctss.df[,3])
     sCTSS <- subset(ctss.df, nTSSs>= expThresh)
     sCTSS.p <- subset(sCTSS, strand>= "+")
     sCTSS.p <- as.matrix(sCTSS.p) #a kludge we'll use for now
     sCTSS.p <- sCTSS.p[complete.cases(sCTSS.p),] 
     my.len <- nrow(sCTSS.p)
     ctss.list.p <- vector(mode="list")
     as.numeric(sCTSS.p[1,1]) -> my.ctss
     j <- 0
     for (i in 1:(my.len-1)) {
            as.numeric(sCTSS.p[i,1]) -> ctss.1 
            as.numeric(sCTSS.p[(i+1),1]) -> ctss.2
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

     sCTSS.m <- subset(sCTSS, strand>= "+")
     sCTSS.m <- as.matrix(sCTSS.m) #a kludge we'll use for now
     sCTSS.m <- sCTSS.m[complete.cases(sCTSS.m),] 
     my.len <- nrow(sCTSS.m)
     ctss.list.m <- vector(mode="list")
     as.numeric(sCTSS.m[1,1]) -> my.ctss
     j <- 0
     for (i in 1:(my.len-1)) {
            as.numeric(sCTSS.m[i,1]) -> ctss.1 
            as.numeric(sCTSS.m[(i+1),1]) -> ctss.2
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
     return(ctss.list)
 }

