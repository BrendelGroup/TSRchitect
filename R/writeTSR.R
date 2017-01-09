#' writeTSR
#' Write TSRs for a specified data set to a tab-delimited file
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param tsrSetType - specifies the set to be written to file. Options are "replicates" or "merged".
#' @param tsrSet - number of the dataset to be processed
#'
#' @return nothing
#' @export

setGeneric(
           name="writeTSR",
           def=function(experimentName, tsrSetType, tsrSet=1) {
               standardGeneric("writeTSR")
    }
    )

setMethod("writeTSR",
          signature(experimentName="tssObject", "character", "numeric"),

          function(experimentName, tsrSetType, tsrSet) {

             message("... writeTSR ...")
             if (tsrSetType=="replicates") {
                 if (tsrSet>length(experimentName@tsrData)) {
                     stop("The value selected for tsrSet exceeds the number of slots in tsrData.")
                 }
                 outfname <- paste("TSRset-", tsrSet, sep="")
                 outfname <- paste(outfname, "txt", sep=".")
                 message("\nThe TSR set for TSS dataset ", tsrSet, " has been written to file ", outfname, "\nin your working directory.")
                 tsr.df <- experimentName@tsrData[[tsrSet]]
             }
             else if (tsrSetType=="merged") {
                 if (length(experimentName@tsrDataMerged)<1) {
                     stop("The @tsrDataMerged slot is currently empty. Please complete the merger before continuing.")
                 }
                 if (tsrSet>length(experimentName@tsrDataMerged)) {
                     stop("The value selected for tsrSet exceeds the number of slots in tsrDataMerged.")
                 }
                 if (tsrSet<length(experimentName@tssCountDataMerged)) {
                     outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                     outfname <- paste(outfname, "txt", sep=".")
                     message("\nThe merged TSR set for TSS dataset ", tsrSet, " has been written to file ", outfname, "\nin your working directory.")
                 }
                 else { # "combined" case
                     outfname <- "TSRsetCombined.txt"
                     message("\nThe combined TSR set derived from all samples has been written to file ", outfname, "\nin your working directory.")
                 }
                 tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
             }
             else {
                 stop("Error: argument tsrSetType to writeTSR() should be either \"replicates\" or \"merged\".")
             }

             write.table(tsr.df, file=outfname, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

             cat("--------------------------------------------------------------------------------\n")
             message(" Done.\n")
          }
          )
