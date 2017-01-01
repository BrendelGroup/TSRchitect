#' findTSR
#' Finds TSRs from a given sequence
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param setToCluster - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param tssSet - number of the dataset to be analyzed
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#'
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="findTSR",
           def=function(experimentName, setToCluster, tssSet=1, tagCountThreshold, clustDist, writeTable=FALSE) {
               standardGeneric("findTSR")
    }
    )

setMethod("findTSR",
          signature(experimentName="tssObject", "character", "numeric", "numeric", "numeric", "logical"),

          function(experimentName, setToCluster, tssSet, tagCountThreshold=1, clustDist, writeTable=FALSE) {

             message("... findTSR ...")
             if (setToCluster=="replicates") {
                 if (tssSet>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds the number of slots in tssCountData.")
                 }
                 tss.df <- experimentName@tssCountData[[tssSet]]
             }
             else if (setToCluster=="merged") {
                 if (length(experimentName@tssCountDataMerged)<1) {
                     stop("The @tssCountDataMerged slot is currently empty. Please complete the merger before continuing.")
                 }
                 if (tssSet>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds the number of slots in tssCountDataMerged.")
                 }
                 tss.df <- experimentName@tssCountDataMerged[[tssSet]]
             }
             else {
                 stop("Error: argument setToCluster to findTSR() should be either \"replicates\" or \"merged\".")
             }

             tsr.list <- tsrCluster(tss.df, minNbrTSSs=tagCountThreshold, minDist=clustDist)
             tsr.DF <- tsrToDF(tsr.list)

             if (writeTable=="TRUE") {
                 if (setToCluster=="replicates") {
                     df.name <- paste("TSRset-", tssSet, sep="")
                     df.name <- paste(df.name, "txt", sep=".")
                     message("\nThe TSR set for TSS dataset ", tssSet, " has been written to file ", df.name, "\nin your working directory.")
                 }
                 else { # setToCluster=="merged" case
                     if (tssSet < length(experimentName@tssCountDataMerged)) {
                         df.name <- paste("TSRsetMerged-", tssSet, sep="")
                         df.name <- paste(df.name, "txt", sep=".")
                         message("\nThe merged TSR set for TSS dataset ", tssSet, " has been written to file ", df.name, "\nin your working directory.")
                     }
                     else {
                         df.name <- "TSRsetCombined.txt"
                         message("\nThe combined TSR set derived from all samples has been written to file ", df.name, "\nin your working directory.")
                     }
                 }
                 write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
             }

	     return(tsr.DF)
             cat("--------------------------------------------------------------------------------\n")
             message(" Done.\n")
          }
          )
