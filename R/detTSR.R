#' detTSR
#' @description An internal function, which is invoked using the user-level
#' function determineTSR that identifies TSRs from the selected tssSet 
#' (Internal function)
#'
#' @param experimentName - a S4 object of class tssObject containing
#' information in slot tssTagData
#' @param tssSetType - specifies the set to be clustered. Options are
#' "replicates" or "merged"
#' @param tssSet - number of the dataset to be analyzed
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#'
#  @keywords internal
#'
#' @return via the user-level function determineTSR, creates a list of
#' GenomicRanges objects containing TSR positions in slot 'tsrData' on
#' the tssObject object
#' 
#' @export
#' @rdname detTSR-methods


setGeneric("detTSR",
    function(experimentName, tssSetType, tssSet=1, tagCountThreshold,
             clustDist)
    standardGeneric("detTSR")
)

#' @rdname detTSR-methods

setMethod("detTSR",
          signature(experimentName="tssObject", "character", "numeric",
                    "numeric", "numeric"),

          function(experimentName, tssSetType, tssSet,
                   tagCountThreshold=1, clustDist) {

              message("... detTSR ...")
              if (tssSetType=="replicates") {
                  if (tssSet>length(experimentName@tssCountData)) {
                      stop("The value selected for tssSet exceeds ",
                           "the number of slots in tssCountData.")
                  }
                  tss.df <- experimentName@tssCountData[[tssSet]]
              }
              else if (tssSetType=="merged") {
                  if (length(experimentName@tssCountDataMerged)<1) {
                      stop("The @tssCountDataMerged slot is currently empty.",
                           "\nPlease complete the merger before continuing.")
                  }
                  if (tssSet>length(experimentName@tssCountData)) {
                      stop("The value selected for tssSet exceeds the number ",
                           "of slots in tssCountDataMerged.")
                  }
                  tss.df <- experimentName@tssCountDataMerged[[tssSet]]
              }
              else {
                  stop("Error: argument tssSetType to detTSR() should be ",
                       "either \"replicates\" or \"merged\".")
              }

              tsr.DF <- tsrCluster(tss.df, minNbrTAGs=tagCountThreshold,
                                   minDist=clustDist)

              message("---------------------------------------------------------\n")
              message(" Done.\n")
              return(tsr.DF)
          }
          )
