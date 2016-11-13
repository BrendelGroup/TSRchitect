#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName - a S4 object of class tssExp containing information in slot tssData
#' @param tssNum - number of the dataset to be analyzed
#' @param nTSSs - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
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

              message("... tsrFind ...")
              if (tssNum>length(expName@expData)) {
                  stop("The value selected for tssNum exceeds the number of slots in tssData.")
              }
              
              tss.mat <- expName@expData[[tssNum]]
              tsr.list <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)
              tsr.DF <- tsrToDF(tsr.list)

              if (writeTable=="TRUE") {
                  df.name <- paste("TSRset-", tssNum, sep="")
                  df.name <- paste(df.name, "txt", sep=".")
                  write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                  message("\nThe TSR set for TSS dataset ", tssNum, " has been written to file ", df.name, "\nin your working directory.")
              }

              expName@tsrData <- tsr.DF
              cat("\n... the TSR dataframe tsrData for dataset ", tssNum, " has been successfully added to\ntssExp object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, envir = parent.frame())              
              message(" Done.\n")
          }
          )
