#' processTSS
#' Finds TSRs from a given sequence
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param parallel - whether to run in parallel or not (logical)
#' @param tssSet - default is "all"; for specific use, specify tssSet number (as character)
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#'
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="processTSS",
           def=function(experimentName, parallel, tssSet, writeTable=FALSE) {
               standardGeneric("processTSS")
    }
    )

setMethod("processTSS",
          signature(experimentName="tssObject", "logical", "character", "logical"),

          function(experimentName, parallel=TRUE, tssSet="all", writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... processTSS ...")
             if (tssSet=="all") {
                 iend <- length(experimentName@replicateIDs)
                 if (parallel==TRUE) {
                     experimentName@tssCountData <- foreach(i=1:iend,.packages="TSRchitect") %dopar% prcTSS(experimentName = experimentName, tssSet = i, writeTable)
                 }
                 else {
                     for (i in 1:iend) {
                         experimentName@tssCountData[[i]] <- prcTSS(experimentName = experimentName, tssSet = i, writeTable)
                     }
                 }
             }
             else {
                 i <- as.numeric(tssSet)
                 if (i>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds the number of slots in tssCountData.")
                 }
                 experimentName@tssCountData[[i]] <- prcTSS(experimentName = experimentName, tssSet = i, writeTable)
             }

             cat("--------------------------------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
