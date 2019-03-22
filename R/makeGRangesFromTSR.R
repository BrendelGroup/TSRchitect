#' @title \strong{makeGRangesFromTSR}
#'
#' @description \code{makeGRangesFromTSR} creates a GRanges object
#' from a specified TSR data set
#'
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@@tssTagData}
#' @param tsrSetType specifies the set to be written to file.
#' Options are "replicates" or "merged". (character)
#' @param tsrSet number of the dataset to be processed (numeric).
#'
#' @return An object of class \emph{GRanges} containing the specified TSR
#' data set. Headers include 'seqnames', 'ranges' (including start and end),
#' 'strand', 'name' (TSR ID) and 'score' (Shape Index/SI) value
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' makeGRangesFromTSR(experimentName=tssObjectExample, tsrSetType="replicates",
#' tsrSet=1)
#'
#' @note For more information on the GRanges class, please visit:
#' http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/GRanges-class.html
#'
#' @export
#' @rdname makeGRangesFromTSR-methods


setGeneric("makeGRangesFromTSR",
    function(experimentName, tsrSetType, tsrSet=1)
    standardGeneric("makeGRangesFromTSR")
)

#' @rdname makeGRangesFromTSR-methods

setMethod("makeGRangesFromTSR",
          signature(experimentName="tssObject", "character", "numeric"),

          function(experimentName, tsrSetType, tsrSet) {

              message("... makeGRangesFromTSR ...")
              if (tsrSetType=="replicates") {
                  if (tsrSet>length(experimentName@tsrData)) {
                      stop("The value selected for tsrSet exceeds the",
                           " number of slots in tsrData.")
                  }
                  tsr.df <- experimentName@tsrData[[tsrSet]]
              }
              else if (tsrSetType=="merged") {
                  if (length(experimentName@tsrDataMerged)<1) {
                      stop("The @tsrDataMerged slot is currently empty.",
                            " Please complete the merger before continuing.")
                  }
                  if (tsrSet>length(experimentName@tsrDataMerged)) {
                      stop("The value selected for tsrSet exceeds the",
                            " number of slots in tsrDataMerged.")
                  }
                  if (tsrSet<length(experimentName@tssCountDataMerged)) {
                      tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
                  }
              else {
                  stop("Error: argument tsrSetType to nakeGRangesFromTSR() ",
                       "should be either \"replicates\" or \"merged\".")
              }
              }

              my.gr <- makeGRangesFromDataFrame
              tsr.df$ID <- paste(tsr.df$seq, tsr.df$start, tsr.df$end,
                                     tsr.df$strand, sep=".")
              bed.df <- tsr.df[, c("seq", "start", "end", "ID",
                                   "tsrSI", "strand")]
              colnames(bed.df) <- c("chrom", "start", "end",
                                    "name", "score", "strand")
              bed.df$start <- as.numeric(as.character(bed.df$start))
              bed.df$end  <- as.numeric(as.character(bed.df$end))
              bed.df$score <- my.score <- as.numeric(as.character(bed.df$score))
              bed.gr <- makeGRangesFromDataFrame(bed.df,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=FALSE,
                                   start.field="start",
                                   end.field="end",
                                   strand.field="strand")
              return(bed.gr)

              message("-------------------------------------------------------\n")
              message(" Done.\n")
          }
          )
