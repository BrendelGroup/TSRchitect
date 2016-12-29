#' tsrFindP
#' Finds TSRs from a given chromosome
#' @param experimentName - a S4 object of class tssObject containing information in slot tssData
#' @param tssSet - number of the dataset to be analyzed
#' @param nTSSs - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param setToCluster - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="tsrFindP",
           def=function(experimentName, tssSet=1, nTSSs, clustDist, setToCluster, writeTable=FALSE) {
               standardGeneric("tsrFindP")
    }
    )

setMethod("tsrFindP",
          signature(experimentName="tssObject", "numeric", "numeric", "numeric", "character", "logical"),

          function(experimentName, tssSet, nTSSs=1, clustDist, setToCluster, writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... tsrFind ...")
             if (setToCluster=="replicates") {
                 if (tssSet>length(experimentName@expData)) {
                     stop("The value selected for tssSet exceeds the number of slots in tssData.")
                 }

                 tss.mat <- experimentName@expData[[tssSet]]
                 tsr.list <- .tsrCluster(tss.mat, expThresh=nTSSs, minDist=clustDist)
                 tsr.DF <- tsrToDF(tsr.list)

                 if (writeTable=="TRUE") {
                     df.name <- paste("TSRset-", tssSet, sep="")
                     df.name <- paste(df.name, "txt", sep=".")
                     write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                     message("\nThe TSR set for TSS dataset ", tssSet, " has been written to file ", df.name, "\nin your working directory.")
                 }

#                experimentName@tsrData[[tssSet]] <- tsr.DF
                 cat("\n... the TSR data frame for dataset ", tssSet, " has been successfully added to\ntssObject object \"", object.name, "\"\n")
		 return(tsr.DF)
              }

              else if (setToCluster=="merged") {
                  if (length(experimentName@expDataMerged)<1) {
                      stop("The @expDataMerged slot is currently empty. Please complete the merger before continuing.")
                  }

                  tsr.list <- vector(mode="list")
                  for (i in 1:length(experimentName@expDataMerged)) {
                      tss.mat <- experimentName@expDataMerged[[i]]
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

                  experimentName@tsrDataMerged <- tsr.list
                  cat("\n... merged TSR data frames have been successfully added to\ntssObject object \"", object.name, "\"\n")
              }
              else {
                  stop("Error: argument setToCluster to tsrFind() should be either \"replicates\" or \"merged\".")
              }
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, envir = parent.frame())
              message(" Done.\n")
          }
          )
