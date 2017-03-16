#' @title \strong{determineTSR}
#' @description \code{determineTSR} Identifies TSRs from
#' entire TSS datasets as specified.
#'
#' @param experimentName an object of class \emph{tssObject}
#' containing information in slot \emph{@@tssTagData}
#' @param n.cores the number of cores to be used for this job.
#' ncores=1 means serial execution of function calls (numeric)
#' @param tsrSetType specifies the set to be clustered.
#' Options are "replicates" or "merged". (character)
#' @param tssSet default is "all"; if a single TSS dataset is desired,
#' specify tssSet number (character)
#' @param tagCountThreshold the number of TSSs required at a given position
#' for it to be considered in TSR identification. (numeric)
#' @param clustDist the maximum distance of TSSs between two
#' TSRs in base pairs. (numeric)
#' @param writeTable specifies whether the output should
#' be written to a table. (logical)
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#'
#' @return creates a list of \linkS4class{GenomicRanges}-containing
#' TSR positions in slot \emph{@@tsrData} of the returned \emph{tssObject}
#' object
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' tssObjectExample <- determineTSR(experimentName=tssObjectExample, n.cores=1,
#' tsrSetType="replicates", tssSet="1", tagCountThreshold=25, clustDist=20,
#' writeTable=FALSE)
#'
#' @note An example similar to this one can be found in the vignette
#' (/inst/doc/TSRchitect.Rmd)

#' @export


setGeneric("determineTSR",
    function(experimentName, n.cores, tsrSetType, tssSet, tagCountThreshold,
             clustDist, writeTable=FALSE)
    standardGeneric("determineTSR")
)

setMethod("determineTSR",
          signature(experimentName="tssObject", "numeric", "character",
                    "character", "numeric", "numeric", "logical"),

          function(experimentName, n.cores=1, tsrSetType=c("replicates",
                   "merged"), tssSet="all", tagCountThreshold=1, clustDist=20,
                   writeTable=FALSE) {

              message("... determineTSR ...")
              fileType <- match.arg(tsrSetType, several.ok=FALSE)
              if (tsrSetType=="replicates") {
                  if (tssSet=="all") {
                      iend <- length(experimentName@tssCountData)
                          multicoreParam <- MulticoreParam(workers=n.cores)
                          FUN  <- function(x) {
                                     detTSR(experimentName=experimentName,
                                     tsrSetType="replicates",
                                     tssSet=x,
                                     tagCountThreshold,
                                     clustDist)
                                 }
                          experimentName@tsrData <- bplapply(1:iend, FUN)
                          if (writeTable==TRUE) {
                              for (i in 1:iend) {
                                   writeTSR(experimentName = experimentName,
                                   tsrSetType="replicates",
                                   tsrSet=i,
                                   fileType="tab")
                               }
                           }
                  }
                  else {
                      i <- as.numeric(tssSet)
                      if (i>length(experimentName@tssCountData)) {
                          stop("The value selected for tssSet",
                               "exceeds the number of slots in tssCountData.")
                      }
                      experimentName@tsrData[[i]] <-
                          detTSR(experimentName = experimentName,
                                 tsrSetType="replicates",
                                 tssSet=i,
                                 tagCountThreshold,
                                 clustDist)
                      if (writeTable==TRUE) {
                          writeTSR(experimentName = experimentName,
                                   tsrSetType="replicates",
                                   tsrSet=i,
                                   fileType="tab")
                      }
                  }
              }

              else if (tsrSetType=="merged") {
                  iend <- length(experimentName@tssCountDataMerged)
                  if (tssSet=="all") {
                      iend <- length(experimentName@tssCountDataMerged)
                          multicoreParam <- MulticoreParam(workers=n.cores)
                          FUN  <- function(x) {
                                     detTSR(experimentName=experimentName,
                                     tsrSetType="merged",
                                     tssSet=x,
                                     tagCountThreshold,
                                     clustDist)
                                     }
                          experimentName@tsrDataMerged <- bplapply(1:iend, FUN)
                          if (writeTable==TRUE) {
                              for (i in 1:iend) {
                                  writeTSR(experimentName =
                                  experimentName,
                                  tsrSetType="merged",
                                  tsrSet=i,
                                  fileType="tab")
                              }
                          }
                  }
                    else {
                          for (i in 1:iend) {
                              experimentName@tsrDataMerged[[i]] <-
                                  detTSR(experimentName = experimentName,
                                         tsrSetType="merged",
                                         tssSet=i,
                                         tagCountThreshold,
                                         clustDist)
                              if (writeTable==TRUE) {
                                  writeTSR(experimentName = experimentName,
                                           tsrSetType="merged",
                                           tsrSet=i,
                                           fileType="tab")
                              }
                          }
                      }
              }
                 else {
                      i <- as.numeric(tssSet)
                      if (i>length(experimentName@tssCountDataMerged)) {
                          stop("The value selected for tssSet exceeds",
                               "the number of slots in tssCountDataMerged.")
                      }
                      experimentName@tsrDataMerged[[i]] <-
                          detTSR(experimentName = experimentName,
                                 tsrSetType="merged",
                                 tssSet=i,
                                 tagCountThreshold,
                                 clustDist)
                      if (writeTable==TRUE) {
                          writeTSR(experimentName = experimentName,
                                   tsrSetType="merged",
                                   tsrSet=i,
                                   fileType="tab")
                      }
                  }
              message("-----------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
