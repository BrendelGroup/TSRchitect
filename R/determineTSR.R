#' @title \strong{determineTSR}
#' @description \code{determineTSR} Identifies TSRs from
#' entire TSS datasets as specified.
#'
#' @param experimentName an object of class \emph{tssObject}
#' containing information in slot \emph{@@tssTagData}
#' @param n.cores the number of cores to be used for this job.
#' ncores=1 means serial execution of function calls (numeric)
#' @param tssSetType specifies the set to be clustered.
#' Options are "replicates" or "merged". (character)
#' @param tssSet default is "all"; if a single TSS dataset is desired,
#' specify tssSet number (character)
#' @param tagCountThreshold the number of TSSs required at a given position
#' for it to be considered in TSR identification. (numeric)
#' @param clustDist the maximum distance of TSSs between two
#' TSRs in base pairs. (numeric)
#' @param writeTable specifies whether the output should
#' be written to a table. (logical)
#' @param mixedorder a logical specifying whether the sequence names should
#' be ordered alphanumerically in the output table ("10" following "9" rather
#' than "1"). (logical)
#'
#' @importFrom BiocParallel bplapply MulticoreParam register
#' @importFrom gtools mixedorder
#'
#' @return creates a list of \linkS4class{GenomicRanges}-containing
#' TSR positions in slot \emph{@@tsrData} of the returned \emph{tssObject}
#' object
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' tssObjectExample <- determineTSR(experimentName=tssObjectExample, n.cores=1,
#' tssSetType="replicates", tssSet="1", tagCountThreshold=25, clustDist=20,
#' writeTable=FALSE)
#'
#' @note An example similar to this one can be found in the vignette
#' (/inst/doc/TSRchitect.Rmd)
#' @export
#' @rdname determineTSR-methods

setGeneric("determineTSR",
    function(experimentName, n.cores, tssSetType, tssSet, tagCountThreshold,
             clustDist, ...)
    standardGeneric("determineTSR")
)

#' @rdname determineTSR-methods

setMethod("determineTSR",
          signature(experimentName="tssObject", "numeric", "character",
                    "character", "numeric", "numeric"),

          function(experimentName, n.cores=1, tssSetType=c("replicates",
                   "merged"), tssSet="all", tagCountThreshold=1, clustDist=20,
                   writeTable=FALSE, mixedorder=FALSE) {

             message("... determineTSR ...")
             if (missing(writeTable)) {writeTable = FALSE}
             if (missing(mixedorder)) {mixedorder = FALSE}
             fileType <- match.arg(tssSetType, several.ok=FALSE)
             if (tssSetType=="replicates") {
                 if (tssSet=="all") {
                     iend <- length(experimentName@tssCountData)
                     fcti <- function(i) {
                                detTSR(experimentName,
                                       tssSetType="replicates",
                                       tssSet=i,
                                       tagCountThreshold,
                                       clustDist)
                             }
                     if (n.cores > 1) {
                         BiocParallel::register(MulticoreParam(workers=n.cores),
							       default=TRUE)
                         experimentName@tsrData <- bplapply(1:iend, fcti)
                     } else {
                         experimentName@tsrData <- lapply(1:iend, fcti)
                     }
                     if (writeTable==TRUE) {
                         for (i in 1:iend) {
                              writeTSR(experimentName,
                                       tsrSetType="replicates",
                                       tsrSet=i,
                                       tsrLabel=paste("TSR",i,sep=""),
                                       mixedorder,
                                       fileType="tab")
                         }
                     }
                 } else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountData)) {
                         stop("The value selected for tssSet",
                              "exceeds the number of slots in tssCountData.")
                     }
                     experimentName@tsrData[[i]] <-
                         detTSR(experimentName = experimentName,
                                tssSetType="replicates",
                                tssSet=i,
                                tagCountThreshold,
                                clustDist)
                     if (writeTable==TRUE) {
                         writeTSR(experimentName,
                                  tsrSetType="replicates",
                                  tsrSet=i,
                                  tsrLabel=paste("TSR",i,sep=""),
                                  mixedorder,
                                  fileType="tab")
                     }
                 }
             } else if (tssSetType=="merged") {
                 iend <- length(experimentName@tssCountDataMerged)
                 if (tssSet=="all") {
                     iend <- length(experimentName@tssCountDataMerged)
                     fcti <- function(i) {
                                detTSR(experimentName=experimentName,
                                       tssSetType="merged",
                                       tssSet=i,
                                       tagCountThreshold,
                                       clustDist)
                             }
                     if (n.cores > 1) {
                         BiocParallel::register(MulticoreParam(workers=n.cores),
							       default=TRUE)
                         experimentName@tsrDataMerged <- bplapply(1:iend, fcti)
                     } else {
                         experimentName@tsrDataMerged <- lapply(1:iend, fcti)
                     }
                     if (writeTable==TRUE) {
                         for (i in 1:iend) {
                              writeTSR(experimentName,
                                       tsrSetType="merged",
                                       tsrSet=i,
                                       tsrLabel=paste("mTSR",i,sep=""),
                                       mixedorder,
                                       fileType="tab")
                         }
                     }
                 } else {
                     i <- as.numeric(tssSet)
                     if (i>length(experimentName@tssCountDataMerged)) {
                         stop("The value selected for tssSet exceeds",
                              "the number of slots in tssCountDataMerged.")
                     }
                     experimentName@tsrDataMerged[[i]] <-
                         detTSR(experimentName = experimentName,
                                tssSetType="merged",
                                tssSet=i,
                                tagCountThreshold,
                                clustDist)
                     if (writeTable==TRUE) {
                         writeTSR(experimentName,
                                  tsrSetType="merged",
                                  tsrSet=i,
                                  tsrLabel=paste("mTSR",i,sep=""),
                                  mixedorder,
                                  fileType="tab")
                     }
                 }
             }
             message("-----------------------------------------------------\n")
             message(" Done.\n")
             return(experimentName)
          }
          )
