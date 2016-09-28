#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName an object of tssExp format containing information in slot tssData
#' @param tssNum the number of the dataset to be analyzed
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssExp object
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom BiocGenerics start end
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @export 

setGeneric(
           name="tsrFind",
           def=function(expName, tssNum, chrName, nTSSs, clustDist) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(expName="tssExp", "numeric", "character", "numeric", "numeric"),

          function(expName, tssNum, c) {

              if (tssNum>length(expName@tssData)) {

                  stop("The value selected exceeds the number of slots in tssData.")

              }

              tss <- clusterChr(tssExp, tssNum, chrName, nTSSs, clustDist)

              message("Creating expression matrix for dataset", tssNum, "...\n\n")
              
              tss.mat <- expressionCTSS(tss)

              message("Clustering TSS expression matrix into TSR regions.\n\n")

              tsr.list <- tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)

              message("Clustering complete.")

              tsr.GR <- tsrToGR(tsr.list, seqname=chrName)

              expName@tsrData <- tsr.GR
                                          
              cat("\nTSRs from", bam.len, "were successfully added to your tssExp object.\n\n")

              assign(object.name, expName, envir = parent.frame())              
              
          }
          )
              


          
