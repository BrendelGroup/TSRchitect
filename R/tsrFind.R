#' tsrFind
#' Finds TSRs from a given chromosome
#' @param expName - a S4 object of class tssExp containing information in slot tssData
#' @param tssNum - number of the dataset to be analyzed
#' @param nTSSs - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param setToCluster - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssExp object
#' @export

setGeneric(
           name="tsrFind",
           def=function(expName, tssNum=1, nTSSs, clustDist, setToCluster, writeTable=FALSE) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(expName="tssExp", "numeric", "numeric", "numeric", "character", "logical"),

          function(expName, tssNum, nTSSs=1, clustDist, setToCluster, writeTable=FALSE) {
             object.name <- deparse(substitute(expName))

             message("... tsrFind ...")
             if (setToCluster=="replicates") {
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

                 expName@tsrData[[tssNum]] <- tsr.DF
                 cat("\n... the TSR data frame for dataset ", tssNum, " has been successfully added to\ntssExp object \"", object.name, "\"\n")
              }

              else if (setToCluster=="merged") {
                  if (length(expName@expDataMerged)<1) {
                      stop("The @expDataMerged slot is currently empty. Please complete the merger before continuing.")
                  }

                  tsr.list <- vector(mode="list")
                  for (i in 1:length(expName@expDataMerged)) {
                      tss.mat <- expName@expDataMerged[[i]]
                      my.tsr <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)
                      tsr.DF <- tsrToDF(my.tsr)
                      tsr.list[[i]] <- tsr.DF

                      if (writeTable=="TRUE") {
                          df.name <- paste("TSRsetMerged-", i, sep="")
                          df.name <- paste(df.name, "txt", sep=".")
                          write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                          message("\nThe merged TSR set for TSS dataset ", i, " has been written to file ", df.name, "\nin your working directory.")
                      }
                  }

                  expName@tsrDataMerged <- tsr.list
                  cat("\n... merged TSR data frames have been successfully added to\ntssExp object \"", object.name, "\"\n")
              }
              else {
                  stop("Error: argument setToCluster to tsrFind() should be either \"replicates\" or \"merged\".")
              }
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, envir = parent.frame())
              message(" Done.\n")
          }
          )
