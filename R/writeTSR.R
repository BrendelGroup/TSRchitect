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
#' @param tsrLabel specifies the label to be used in the name column of
#' BED format output
#' @param mixedorder a logical specifying whether the sequence names should
#' be ordered alphanumerically ("10" following "9" rather than "1"). (logical)
#' @param fileType the format of the file to be written.
#' Possible choices are "tab" for tab-delimited output, "bed" for BED format,
#' and "gff" for for GFF3 format (character).
#'
#' @return A table containing the specified TSR data set that
#' is to be written to your working directory.
#'
#' @import     BiocGenerics
#' @importFrom gtools mixedorder
#' @importFrom stats median
#' @importFrom rtracklayer export.bed
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' writeTSR(experimentName=tssObjectExample, tsrSetType="replicates",
#'          tsrSet=1, tsrLabel="TSRsample1_", mixedorder=FALSE, fileType="tab")
#'
#' @note The .bed file written adheres to the standard six-column BED format,
#' while "tab" format is identical to that of the data.frames containing TSR
#' data.
#' @note For more information on the BED format, please visit
#' https://genome.ucsc.edu/FAQ/FAQformat#format1
#'
#' @export
#' @rdname writeTSR-methods


setGeneric("writeTSR",
    function(experimentName, tsrSetType, tsrSet=1, tsrLabel="TSR_",
             mixedorder=FALSE,fileType="tab")
    standardGeneric("writeTSR")
)

#' @rdname writeTSR-methods

setMethod("writeTSR",
          signature(experimentName="tssObject", "character", "numeric",
                    "character", "logical", "character"),

          function(experimentName, tsrSetType, tsrSet,
                   tsrLabel, mixedorder, fileType) {

              message("... writeTSR ...")
              if (tsrSetType=="replicates") {
                  if (tsrSet>length(experimentName@tsrData)) {
                      stop("The value selected for tsrSet exceeds the",
                           " number of slots in tsrData.")
                  }
                  outfname <- paste("TSRset-", tsrSet, sep="")
                  if (fileType == "tab") {
                      outfname <- paste(outfname, "tab", sep=".")
                  } else if (fileType == "bed") {
                      outfname <- paste(outfname, "bed", sep=".")
                  } else if (fileType == "gff") {
                      outfname <- paste(outfname, "gff", sep=".")
                  } else {
                      stop("Unknown fileType selected for writeTSR.",
                           " Please check.")
                  }
                  message("\nThe TSR set for TSS dataset ", tsrSet,
                          " will be written to file ",
                          outfname, "\nin your working directory.")
                  if (!missing(mixedorder) & mixedorder == TRUE) {
                    tsr.df <- experimentName@tsrData[[tsrSet]][mixedorder(
                                   experimentName@tsrData[[tsrSet]]$seq),]
                  } else {
                    tsr.df <- experimentName@tsrData[[tsrSet]]
                  }
              } else if (tsrSetType=="merged") {
                  if (length(experimentName@tsrDataMerged)<1) {
                      stop("The @tsrDataMerged slot is currently empty.",
                           " Please complete the merger before continuing.")
                  }
                  if (tsrSet>length(experimentName@tsrDataMerged)) {
                      stop("The value selected for tsrSet exceeds the",
                           " number of slots in tsrDataMerged.")
                  }
                  if (tsrSet<length(experimentName@tssCountDataMerged)) {
                      outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                      if (fileType == "tab") {
                          outfname <- paste(outfname, "tab", sep=".")
                      } else if (fileType == "bed") {
                          outfname <- paste(outfname, "bed", sep=".")
                      } else if (fileType == "gff") {
                          outfname <- paste(outfname, "gff", sep=".")
                      } else {
                          stop("Unknown fileType selected for writeTSR.",
                               " Please check.")
                      }
                      message("\nThe merged TSR set for TSS dataset ", tsrSet,
                      " will be written to file ", outfname,
                      "\nin your working directory.")
                  } else { # "combined" case
                      if (fileType == "tab") {
                          outfname <- "TSRsetCombined.tab"
                      } else if (fileType == "bed") {
                          outfname <- "TSRsetCombined.bed"
                      } else if (fileType == "gff") {
                          outfname <- "TSRsetCombined.gff"
                      } else {
                          stop("Unknown fileType selected for writeTSR.",
                               " Please check.")
                      }
                      message("\nThe combined TSR set derived from all samples",
                              " will be written to file ", outfname,
                              "\nin your working directory.")
                  }
                  if (!missing(mixedorder) & mixedorder == TRUE) {
                    tsr.df <- experimentName@tsrDataMerged[[tsrSet]][mixedorder(
                                   experimentName@tsrDataMerged[[tsrSet]]$seq),]
                  } else {
                    tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
                  }
              } else {
                  stop("Error: argument tsrSetType to writeTSR() should be",
                       " either \"replicates\" or \"merged\".")
              }

              if (fileType == "tab") {
                  write.table(format(tsr.df,scientific=FALSE), file=outfname,
                              col.names=FALSE, row.names=FALSE, sep="\t",
                              quote=FALSE)
              } else if (fileType == "bed") {
                  tsr.df$ID <- paste(tsrLabel,which(tsr.df$seq != ""),sep="_")
                  m <- pmax(median(tsr.df$nTAGs),1)
                  tsr.df$score <- round(100*pmin(tsr.df$nTAGs/m,10),0)
                  out.df <- tsr.df[, c("seq", "start", "end", "ID",
                                       "score", "strand")]
                  colnames(out.df) <- c("chrom", "start", "end","name",
                                        "score", "strand")
                  export.bed(out.df,con=outfname)
              } else { # fileType == "gff") 
                  tsr.df$source <- rep("TSRchitect",nrow(tsr.df))
                  tsr.df$type <- rep("TSR",nrow(tsr.df))
                  m <- pmax(median(tsr.df$nTAGs),1)
                  tsr.df$score <- round(100*pmin(tsr.df$nTAGs/m,10),0)
                  tsr.df$phase <- rep(".",nrow(tsr.df))
                  tsr.df$ID <- paste("ID=",tsrLabel,"_",
                                     which(tsr.df$seq != ""),";",sep="")
                  out.df <- tsr.df[with(tsr.df,order(tsr.df$seq,tsr.df$start,
                                                     tsr.df$strand)),
                                   c("seq", "source", "type", "start", "end",
                                     "score", "strand", "phase", "ID")]
                  write.table(format(out.df,scientific=FALSE,trim=TRUE),
                              file=outfname, col.names=FALSE, row.names=FALSE,
                              sep="\t", quote=FALSE)
              }

              message("---------------------------------------------------------\n")
              message(" Done.\n")
          }
          )
