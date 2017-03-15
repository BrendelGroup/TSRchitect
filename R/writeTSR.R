#' @title \strong{writeTSR}
#'
#' @description \code{writeTSR} writes identified TSRs
#' from a specified data set to a file in either tab or BED formats
#'
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@@tssTagData}
#' @param tsrSetType specifies the set to be written to file.
#' Options are "replicates" or "merged". (character)
#' @param tsrSet number of the dataset to be processed (numeric).
#' @param fileType the format of the file to be written.
#' Possible choices are "tab" for tab-delimited output and
#' "bed" for BED format (character).
#'
#' @return A table containing the specified TSR data set that
#' is to be written to your working directory.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer export.bed
#' @importFrom utils write.table
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' writeTSR(experimentName=tssObjectExample, tsrSetType="replicates",
#' tsrSet=1, fileType="tab")
#'
#' @note The .bed file written adheres to the standard six-column BED format,
#' while "tab" format is identical to that of the data.frames containing TSR
#' data.
#' @note For more information on the BED format, please visit
#' https://genome.ucsc.edu/FAQ/FAQformat#format1
#'
#' @export


setGeneric("writeTSR",
    function(experimentName, tsrSetType, tsrSet=1, fileType="tab")
    standardGeneric("writeTSR")
)

setMethod("writeTSR",
          signature(experimentName="tssObject", "character", "numeric",
		    "character"),

          function(experimentName, tsrSetType, tsrSet, fileType) {

              message("... writeTSR ...")
              if (tsrSetType=="replicates") {
                  if (tsrSet>length(experimentName@tsrData)) {
                      stop("The value selected for tsrSet exceeds the",
                           "number of slots in tsrData.")
                  }
                  outfname <- paste("TSRset-", tsrSet, sep="")
                  if (fileType == "tab") {
                      outfname <- paste(outfname, "tab", sep=".")
                  }
                  else if (fileType == "bed") {
                      outfname <- paste(outfname, "bed", sep=".")
                  }
                  else {
                      stop("Unknown fileType selected for writeTSR.",
                           "Please check.")
                  }
                  message("\nThe TSR set for TSS dataset ", tsrSet,
                          " has been written to file ",
                         outfname, "\nin your working directory.")
                  tsr.df <- experimentName@tsrData[[tsrSet]]
              }
              else if (tsrSetType=="merged") {
                  if (length(experimentName@tsrDataMerged)<1) {
                      stop("The @tsrDataMerged slot is currently empty.",
                            "Please complete the merger before continuing.")
                  }
                  if (tsrSet>length(experimentName@tsrDataMerged)) {
                      stop("The value selected for tsrSet exceeds the",
                            "number of slots in tsrDataMerged.")
                  }
                  if (tsrSet<length(experimentName@tssCountDataMerged)) {
                      outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                      if (fileType == "tab") {
                          outfname <- paste(outfname, "tab", sep=".")
                      }
                      else if (fileType == "bed") {
                          outfname <- paste(outfname, "bed", sep=".")
                      }
                      else {
                          stop("Unknown fileType selected for writeTSR.",
                               "Please check.")
                      }
                      message("\nThe merged TSR set for TSS dataset ", tsrSet,
                      " has been written to file ", outfname,
                      "\nin your working directory.")
                  }
                  else { # "combined" case
                      if (fileType == "tab") {
                          outfname <- "TSRsetCombined.tab"
                      }
                      else if (fileType == "bed") {
                          outfname <- "TSRsetCombined.bed"
                      }
                      else {
                          stop("Unknown fileType selected for writeTSR.
                          Please check.")
                      }
                      message("\nThe combined TSR set derived from all samples",
                              " has been written to file ", outfname,
                              "\nin your working directory.")
                  }
                  tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
              }
              else {
                  stop("Error: argument tsrSetType to writeTSR() should be",
                       "either \"replicates\" or \"merged\".")
              }

              if (fileType == "tab") {
                  write.table(tsr.df, file=outfname, col.names=FALSE,
                  row.names=FALSE, sep="\t", quote=FALSE)
              }
              else {
                  tsr.df$ID <- paste(tsr.df$seq, tsr.df$start, tsr.df$end,
                               tsr.df$strand, sep=".")
                  bed.df <- tsr.df[, c("seq", "start", "end", "ID",
                             "shapeIndex", "strand")]
                  colnames(bed.df) <- c("chrom", "start", "end",
                                      "name", "score", "strand")
                  bed.df$start <- bed.df$start
                  bed.gr <- makeGRangesFromDataFrame(bed.df,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           start.field="start",
                                           end.field="end",
                                           strand.field="strand")
                  export.bed(bed.gr, con=outfname)
              }

              message("---------------------------------------------------------\n")
              message(" Done.\n")
          }
          )
