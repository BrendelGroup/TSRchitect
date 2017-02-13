#' @title \strong{processTSS}
#' @description \code{processTSS} calulates the number of observed reads
#' at a given TSS coordinate across an entire dataset.
#' 
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@tssTagData}
#' @param parallel if TRUE, the job will be run in parallel (logical)
#' @param n.cores the number of cores to be used for this job.
#' Ignored if 'parallel=FALSE' (numeric)
#' @param tssSet default is "all"; to select a single \emph{tssSet},
#' specify it (as character)
#' @param writeTable specifies whether the output should be written
#' to a file. (logical)
#'
#' @importFrom BiocParallel bplapply MulticoreParam
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
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' 
#' @export

setGeneric(
           name="processTSS",
    def=function(experimentName, parallel, n.cores, tssSet, writeTable) {
               standardGeneric("processTSS")
    }
    )

setMethod("processTSS",
          signature("tssObject", "logical", "numeric", "character", "logical"),

          function(experimentName, parallel=TRUE, n.cores=1,
                   tssSet="all", writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... processTSS ...")
             if (tssSet=="all") {
                 iend <- length(experimentName@tssTagData)
                 if (parallel==TRUE) {
                     multicoreParam <- MulticoreParam(workers=n.cores)
                     FUN  <- function(x) {
                            prcTSS(experimentName, tssSet=x,
                            writeTable=writeTable)
                            }
                     experimentName@tssCountData <- bplapply(1:iend, FUN)
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
                 if (i > length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds",
                          "the number of slots in tssCountData.")
                 }
                 experimentName@tssCountData[[i]] <- prcTSS(experimentName =
                                                            experimentName,
                                                            tssSet = i,
                                                            writeTable)
             }
             cat("--------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
