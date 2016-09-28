#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName an object of tssExp format containing information in slot tssData
#' @param tssNum the number of the dataset to be analyzed
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssExp object
#' @export 

setGeneric(
           name="tsrFind",
           def=function(expName, tssNum, chrName, nTSSs, clustDist) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(expName="tssExp", "numeric", "character", "numeric", "numeric"),

          function(expName, tssNum, chrName, nTSSs, clustDist) {

              object.name <- deparse(substitute(expName))
              
              message("\nInitiated TSR finding.")

              if (tssNum>length(expName@tssData)) {

                  stop("The value selected exceeds the number of slots in tssData.")

              }

              tss <- tssChr(expName, tssNum, chrName)

              message("Creating expression matrix for dataset ", tssNum, "...\n")
              
              tss.mat <- expressionCTSS(tss)

              message("Clustering TSS expression matrix into TSR regions.\n")

              tsr.list <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)

              message("Clustering complete.")

              tsr.GR <- tsrToGR(tsr.list, seqname=chrName)

              expName@tsrData <- tsr.GR
                                          
              cat("\nTSRs from were successfully added to your tssExp object.\n")

              assign(object.name, expName, envir = parent.frame())              
              
          }
          )
              


          
