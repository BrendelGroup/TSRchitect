#' Partitions, then clusters tss data by chromosome.
#' @param expName an object of class tssExp containing the tss data
#' @param tssNum the slot number of the tss data to be clustered
#' @param chrName the name of the chromosome to select
#' @return a list object containing xmeans-clustered TSS data for a single chromosome (plus and minus strands)
#' @importFrom BiocGenerics strand start
#' @importFrom GenomeInfoDb seqnames
#' @export

setGeneric(
           name="clusterChr",
           def=function(expName, tssNum, chrName) {
               standardGeneric("clusterChr")
    }
    )

setMethod("clusterChr",
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

                  tss.plus.clustered <- xmeans(tss.plus.vec, ik=100)

                  chr.list$plus <- tss.plus.clustered

                  chr.list$plus$tss <- tss.plus.vec

                  # now for the TSSs on the minus strand

                  tss.minus <- tss.total[strand(tss.total)=="-",]

                  tss.minus.vec <- start(tss.minus)

                  tss.minus.clustered <- xmeans(tss.minus.vec, ik=100, iter.max=10)

                  chr.list$minus <- tss.minus.clustered

                  chr.list$minus$tss <- tss.minus.vec

                  tss.total <- NULL

                  return(chr.list)
          }
           )   

