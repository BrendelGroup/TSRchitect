#' @title \strong{processTSS}
#' @description \code{processTSS} calulates the number of observed reads
#' at a given TSS coordinate across an entire dataset.
#'
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@tssTagData}
#' @param n.cores the number of cores to be used for this job.
#' n.cores=1 means serial execution of function calls (numeric)
#' @param tssSet default is "all"; to select a single \emph{tssSet},
#' specify it (as character)
#' @param writeTable specifies whether the output should be written
#' to a file. (logical)
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#'
#' @return Creates a list of \linkS4class{GenomicRanges} containing TSS
#' positions in slot \emph{tssTagData} of the returned \emph{tssObject}.
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' tssObjectExample <- processTSS(experimentName=tssObjectExample, n.cores=1,
#' tssSet="all", writeTable=FALSE)
#'
#' @note Note that the \emph{tssSet} parameter must be of class
#' \emph{character}, even when selecting an individual dataset.
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#'
#' @export
#' @rdname processTSS-methods


setGeneric("processTSS",
    function(experimentName, n.cores, tssSet, writeTable)
    standardGeneric("processTSS")
)

#' @rdname processTSS-methods

setMethod("processTSS",
          signature("tssObject", "numeric", "character", "logical"),

          function(experimentName, n.cores=1, tssSet="all", writeTable=FALSE) {

             message("... processTSS ...")
             if (tssSet=="all") {
                 iend <- length(experimentName@tssTagData)
                 funi <- function(i) {
                            prcTSS(experimentName, n.cores, tssSet=i,
                                   writeTable)
                         }
                 if (n.cores > 1) {
		     BiocParallel::register(MulticoreParam(workers=n.cores),
							   default=TRUE)
                     experimentName@tssCountData <- bplapply(1:iend, funi)
                 }
		 else {
                     experimentName@tssCountData <- lapply(1:iend, funi)
                 }
             }
             else {
                 i <- as.numeric(tssSet)
                 if (i > length(experimentName@tssTagData)) {
                     stop("The value selected for tssSet exceeds ",
                          "the number of slots in tssTagData.")
                 }
                 experimentName@tssCountData[[i]] <-
		     prcTSS(experimentName, n.cores=1, tssSet = i, writeTable)
              }
              message("-----------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
