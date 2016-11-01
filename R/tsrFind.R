#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName an object of tssExp format containing information in slot tssData
#' @param tssNum the number of the dataset to be analyzed
#' @param nTSSs the number of TSSs required at a given position
#' @param clustDist the maximum distance of TSSs between two TSRs (in base pairs)
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssExp object
#' @export 

setGeneric(
           name="tsrFind",
           def=function(expName, tssNum, nTSSs, clustDist) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(expName="tssExp", "numeric", "numeric", "numeric"),

          function(expName, tssNum, nTSSs, clustDist) {
              object.name <- deparse(substitute(expName))
              message("\nInitiated TSR finding.")

              if (tssNum>length(expName@expData)) {
                  stop("The value selected exceeds the number of slots in tssData.")
              }
              
              tss.mat <- expName@expData[[tssNum]]
              message("Clustering TSS expression matrix into TSR regions.\n")
              tsr.list <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)
              message("Clustering complete.\n")
              tsr.DF <- tsrToDF(tsr.list)
              expName@tsrData <- tsr.DF
              message("\nTSRs from were successfully added to your tssExp object.\n")
              assign(object.name, expName, envir = parent.frame())              
          }
          )
              


          
