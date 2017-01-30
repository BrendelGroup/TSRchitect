#' @title \strong{processTSS}
#' @description \code{processTSS} calulates the number of observed reads
#' at a given TSS coordinate across an entire dataset.
#' 
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@tssTagData}
#' @param parallel whether to run in parallel or not (logical)
#' @param tssSet default is "all"; for specific use, specify \emph{tssSet}
#' number (as character)
#' @param writeTable specifies whether the output should be written
#' to a table. (logical)
#' 
#' @return Creates a list of \linkS4class{GenomicRanges} containing TSS
#' positions in slot \emph{tssTagData} on the \emph{tssObject}.
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' processTSS(experimentName=tssObjectExample, parallel=FALSE, tssSet="all",
#' writeTable=FALSE)
#' 
#' @note Note that the \emph{tssSet} parameter must be of class
#' \emph{character}, even when selecting an individual dataset.
#' @note An example similar to the one provided can be found in
#' \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd).
#' 
#' @export

setGeneric(
           name="processTSS",
           def=function(experimentName, parallel, tssSet, writeTable=FALSE) {
               standardGeneric("processTSS")
    }
    )

setMethod("processTSS",
          signature(experimentName="tssObject", "logical", "character",
                    "logical"),

          function(experimentName, parallel=TRUE, tssSet="all",
                   writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... processTSS ...")
             if (tssSet=="all") {
                 iend <- length(experimentName@replicateIDs)
                 if (parallel==TRUE) {
                     experimentName@tssCountData <- foreach(i=1:iend,
                     .packages="TSRchitect") %dopar% 
                     prcTSS(experimentName = experimentName,
                                     tssSet = i, writeTable)
                 }
                 else {
                     for (i in 1:iend) {
 experimentName@tssCountData[[i]] <- prcTSS(experimentName=experimentName,
                                            tssSet = i, writeTable=TRUE)
                     }
                 }
             }
             else {
                 i <- as.numeric(tssSet)
                 if (i>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds",
                          "the number of slots in tssCountData.")
                 }
 experimentName@tssCountData[[i]] <- prcTSS(experimentName = experimentName,
                                      tssSet = i, writeTable)
             }

             cat("--------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
