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
           def=function(expName, tssNum, nTSSs, clustDist, writeTable) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(expName="tssExp", "numeric", "numeric", "numeric", "logical"),

          function(expName, tssNum, nTSSs, clustDist, writeTable) {
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

              if (writeTable=="TRUE") {
              df.name <- paste("CTSS", tssNum, sep="")
              df.name <- paste(df.name, "txt", sep=".")
              write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
              message("\nA data frame containing TSRs was written to your current working directory.")
              }

              expName@tsrData <- tsr.DF
              
              message("\nTSRs were successfully added to your tssExp object.\n")
              assign(object.name, expName, envir = parent.frame())              
          }
          )
              


          
