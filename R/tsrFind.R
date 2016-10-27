#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName an object of tssExp format containing information in slot tssData
#' @param tssNum the number of the dataset to be analyzed
#' @param clustDist the maximum distance of TSSs between two TSRs (in base pairs)
#' @param nTSSs the number of TSSs required at a given position
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

              if (tssNum>length(expName@tssData)) {
                  stop("The value selected exceeds the number of slots in tssData.")
              }
              
              tss <- acquireTSS(expName, tssNum)

              message("\nCreating expression matrix for dataset ", tssNum, "...\n")

              df.name <- paste("CTSS", tssNum, sep="")
              df.name <- paste(df.name, "txt", sep=".")
              tss.mat <- expressionCTSS(tss, writeDF=TRUE, dfName=df.name)

              message("Clustering TSS expression matrix into TSR regions.\n")

              tsr.list <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)

              message("Clustering complete.")

              tsr.GR <- tsrToGR(tsr.list)
              expName@tsrData <- tsr.GR
              
              message("\nTSRs from were successfully added to your tssExp object.\n")
              
              assign(object.name, expName, envir = parent.frame())              
          }
          )
              


          
