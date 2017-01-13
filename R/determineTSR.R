#' determineTSR
#' Finds TSRs from a given sequence
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param paralel - whether to run in parallel or not (logical)
#' @param tsrSetType - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param tssSet - default is "all"; for specific use, specify tssSet number (as character)
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#'
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="determineTSR",
           def=function(experimentName, parallel, tsrSetType, tssSet, tagCountThreshold, clustDist, writeTable=FALSE) {
               standardGeneric("determineTSR")
    }
    )

setMethod("determineTSR",
          signature(experimentName="tssObject", "logical", "character", "character", "numeric", "numeric", "logical"),

          function(experimentName, parallel=TRUE, tsrSetType, tssSet="all", tagCountThreshold=1, clustDist=20, writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... determineTSR ...")
             if (tsrSetType=="replicates") {
                 if (tssSet=="all") {
                     iend <- length(experimentName@tssCountData)
                     if (parallel==TRUE) {
                         experimentName@tsrData <- foreach(i=1:iend,.packages="TSRchitect") %dopar% detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet = i, tagCountThreshold, clustDist)
                         if (writeTable=="TRUE") {
                             foreach(i=1:iend,.packages="TSRchitect") %dopar% writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet = i)
                         }
                     }
                     else {
                         for (i in 1:iend) {
                             experimentName@tsrData[[i]] <- detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet = i, tagCountThreshold, clustDist)
                             if (writeTable=="TRUE") {
                                 writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet=i, filetype="tab")
                             }
                         }
                     }
                 }
                 else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountData)) {
                         stop("The value selected for tssSet exceeds the number of slots in tssCountData.")
                     }
                     experimentName@tsrData[[i]] <- detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet = i, tagCountThreshold, clustDist)
                     if (writeTable=="TRUE") {
                         writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet=i, filetype="tab")
                     }
                 }
             }

             else if (tsrSetType=="merged") {
                 if (tssSet=="all") {
                     iend <- length(experimentName@tssCountDataMerged)
                     if (parallel==TRUE) {
                         experimentName@tsrDataMerged <- foreach(i=1:iend,.packages="TSRchitect") %dopar% detTSR(experimentName = experimentName, tsrSetType="merged", tssSet = i, tagCountThreshold, clustDist)
                         if (writeTable=="TRUE") {
                             foreach(i=1:iend,.packages="TSRchitect") %dopar% writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, filetype="tab")
                         }
                     }
                     else {
                         for (i in 1:iend) {
                             experimentName@tsrDataMerged[[i]] <- detTSR(experimentName = experimentName, tsrSetType="merged", tssSet = i, tagCountThreshold, clustDist)
                             if (writeTable=="TRUE") {
                                 writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, filetype="tab")
                             }
                         }
                     }
                 }
                 else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountDataMerged)) {
                         stop("The value selected for tssSet exceeds the number of slots in tssCountDataMerged.")
                     }
                     experimentName@tsrDataMerged[[i]] <- detTSR(experimentName = experimentName, tsrSetType="merged", tssSet = i, tagCountThreshold, clustDist)
                     if (writeTable=="TRUE") {
                         writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, filetype="tab")
                     }
                 }
             }

             else {
                 stop("Error: argument tsrSetType to determineTSR() should be either \"replicates\" or \"merged\".")
             }
             cat("--------------------------------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
