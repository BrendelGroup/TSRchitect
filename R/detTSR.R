#' detTSR
#' Finds TSRs from a given sequence
#'
#' @param experimentName - a S4 object of class tssObject containing
#' information in slot tssTagData
#' @param tsrSetType - specifies the set to be clustered.
#' Options are "replicates" or "merged".
#' @param tssSet - number of the dataset to be analyzed
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#'
#' @keywords internal
#'
#' @return creates a list of GenomicRanges containing TSR positions in slot
#' 'tsrData' on your tssObject object

setGeneric(
           name="detTSR",
    def=function(experimentName, tsrSetType, tssSet=1,
        tagCountThreshold, clustDist) {
               standardGeneric("detTSR")
    }
    )

setMethod("detTSR",
          signature(experimentName="tssObject", "character", "numeric",
                    "numeric", "numeric"),

          function(experimentName, tsrSetType, tssSet,
                   tagCountThreshold=1, clustDist) {

             message("... detTSR ...")
             if (tsrSetType=="replicates") {
                 if (tssSet>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds /
 the number of slots in tssCountData.")
                 }
                 tss.df <- experimentName@tssCountData[[tssSet]]
             }
             else if (tsrSetType=="merged") {
                 if (length(experimentName@tssCountDataMerged)<1) {
                     stop("The @tssCountDataMerged slot is currently empty. /
 Please complete the merger before continuing.")
                 }
                 if (tssSet>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds the number /
 of slots in tssCountDataMerged.")
                 }
                 tss.df <- experimentName@tssCountDataMerged[[tssSet]]
             }
             else {
                 stop("Error: argument tsrSetType to detTSR() should be /
 either \"replicates\" or \"merged\".")
             }

             tsr.list <- tsrCluster(tss.df, minNbrTSSs=tagCountThreshold,
                                    minDist=clustDist)
             tsr.DF <- tsrToDF(tsr.list)

             return(tsr.DF)
             cat("---------------------------------------------------------\n")
             message(" Done.\n")
          }
          )
