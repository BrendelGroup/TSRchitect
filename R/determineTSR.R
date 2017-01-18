#' @title \emph{determineTSR}
#' @description \code{determineTSR} Identifies TSRs from entire TSS datasets as specified
#' @param experimentName an object of class \emph{tssObject} containing information in slot \emph{@@tssTagData}
#' @param parallel if TRUE, the analysis is run in parallel (logical)
#' @param tsrSetType specifies the set to be clustered. Options are "replicates" or "merged". (character)
#' @param tssSet default is "all"; if a single TSS dataset is desired, specify tssSet number (character)
#' @param tagCountThreshold the number of TSSs required at a given position for it to be considered in TSR identification. (numeric)
#' @param clustDist the maximum distance of TSSs between two TSRs in base pairs. (numeric)
#' @param writeTable specifies whether the output should be written to a table. (logical)
#' @examples 
#' @note
#' @return creates a list of \linkS4class{GenomicRanges}-containing TSR positions in slot \emph{@@tsrData} on the \emph{tssObject} object
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
                         experimentName@tsrData <- foreach(i=1:iend,.packages="TSRchitect") %dopar% detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet=i, tagCountThreshold, clustDist)
                         if (writeTable=="TRUE") {
                             foreach(i=1:iend,.packages="TSRchitect") %dopar% writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet=i, fileType="tab")
                         }
                     }
                     else {
                         for (i in 1:iend) {
                             experimentName@tsrData[[i]] <- detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet=i, tagCountThreshold, clustDist)
                             if (writeTable=="TRUE") {
                                 writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet=i, fileType="tab")
                             }
                         }
                     }
                 }
                 else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountData)) {
                         stop("The value selected for tssSet exceeds the number of slots in tssCountData.")
                     }
                     experimentName@tsrData[[i]] <- detTSR(experimentName = experimentName, tsrSetType="replicates", tssSet=i, tagCountThreshold, clustDist)
                     if (writeTable=="TRUE") {
                         writeTSR(experimentName = experimentName, tsrSetType="replicates", tsrSet=i, fileType="tab")
                     }
                 }
             }

             else if (tsrSetType=="merged") {
                 if (tssSet=="all") {
                     iend <- length(experimentName@tssCountDataMerged)
                     if (parallel==TRUE) {
                         experimentName@tsrDataMerged <- foreach(i=1:iend,.packages="TSRchitect") %dopar% detTSR(experimentName = experimentName, tsrSetType="merged", tssSet=i, tagCountThreshold, clustDist)
                         if (writeTable=="TRUE") {
                             foreach(i=1:iend,.packages="TSRchitect") %dopar% writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, fileType="tab")
                         }
                     }
                     else {
                         for (i in 1:iend) {
                             experimentName@tsrDataMerged[[i]] <- detTSR(experimentName = experimentName, tsrSetType="merged", tssSet=i, tagCountThreshold, clustDist)
                             if (writeTable=="TRUE") {
                                 writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, fileType="tab")
                             }
                         }
                     }
                 }
                 else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountDataMerged)) {
                         stop("The value selected for tssSet exceeds the number of slots in tssCountDataMerged.")
                     }
                     experimentName@tsrDataMerged[[i]] <- detTSR(experimentName = experimentName, tsrSetType="merged", tssSet=i, tagCountThreshold, clustDist)
                     if (writeTable=="TRUE") {
                         writeTSR(experimentName = experimentName, tsrSetType="merged", tsrSet=i, fileType="tab")
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
