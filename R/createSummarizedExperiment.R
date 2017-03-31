#' @title \strong{createSummarizedExperiment}
#'
#' @description \code{createSummarizedExperiment} creates and reutrns
#' a \emph{SummarizedExperiment} object from the tag counts on a selected 
#' TSR dataset.
#'
#' @param experimentName an S4 object of class \emph{tssObject} containing
#' information in slot \emph{@@tssTagData}
#' @param tsrSetType specifies the TSR set to be converted into a
#' \emph{SummarizedExperiment} object. Options are "replicates" or "merged".
#' (character)
#' @param tsrSet number of the dataset to be processed (numeric).
#' @param samplePrefix the prefix (or prefixes) that match the sample
#' identifiers in the tsrData column. (character)
#'
#' @return A table containing the specified TSR data set that
#' is to be written to your working directory.
#'
#' @import     BiocGenerics
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' createSumarizedExperiment(tssObjectExample, tsrSetType="merged", tsrSet=1,
#' samplePrefix=c("sample1","sample2"))
#'
#' @note For more information on the SummarizedExperiment class, please visit
#' https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
#'
#' @export


setGeneric("createSummarizedExperiment",
    function(experimentName, tsrSetType="merged", tsrSet=1, samplePrefix)
    standardGeneric("createSummarizedExperiment")
)

setMethod("createSummarizedExperiment",
          signature("tssObject", "character", "numeric",
                    "character"),

          function(experimentName, tsrSetType, tsrSet, samplePrefix) {

              message("... createSummarizedExperiment ...")

              if (tsrSetType=="replicates") {
                  if (tsrSet>length(experimentName@tsrData)) {
                      stop("The value selected for tsrSet exceeds the",
                           " number of slots in tsrData.")
                  }
                  my.tsrs <- experimentName@tsrData[[tsrSet]]
                  if (length(samplePrefix) == 1) {
                      my.names <- names(my.tsrs)
                      this.ind <- grep(pattern=samplePrefix, x=my.names)
                      if (length(this.ind)<1) {
                          stop("The samplePrefix value you entered doesn't",
                               " match any columns. Please check.")
                      }
                      my.counts <- my.tsrs[,this.ind]
                      my.ma <- as.matrix(my.counts)
                  }
                  if (length(samplePrefix) > 1) {
                      my.names <- names(my.tsrs)
                      pre.len <- length(samplePrefix)
                      this.ind <- vector(mode="integer")
                      for (i in 1:pre.len) {
                          this.pre <- samplePrefix[i]
                          my.ind <- grep(pattern=this.pre, x=my.names)
                          this.ind <- c(this.ind, my.ind)
                          }
                          if (length(this.ind)<1) {
                              stop("The samplePrefix value you entered doesn't",
                               " match any columns. Please check.")
                          }
                          my.counts <- my.tsrs[,this.ind]
                          my.ma <- as.matrix(my.counts)
                  }
              }

              if (tsrSetType=="merged") {
                  if (tsrSet>length(experimentName@tsrDataMerged)) {
                      stop("The value selected for tsrSet exceeds the",
                           " number of slots in tsrData.")
                  }
                  my.tsrs <- experimentName@tsrDataMerged[[tsrSet]]
                  if (length(samplePrefix) == 1) {
                      my.names <- names(my.tsrs)
                      this.ind <- grep(pattern=samplePrefix, x=my.names)
                      if (length(this.ind)<1) {
                         stop("The samplePrefix value you entered doesn't",
                             " match any columns. Please check.")
                      }
                      my.counts <- my.tsrs[,this.ind]
                      my.ma <- as.matrix(my.counts)
                  }
                  if (length(samplePrefix) > 1) {
                      my.names <- names(my.tsrs)
                      pre.len <- length(samplePrefix)
                      this.ind <- vector(mode="integer")
                      for (i in 1:pre.len) {
                          this.pre <- samplePrefix[i]
                          my.ind <- grep(pattern=this.pre, x=my.names)
                          this.ind <- c(this.ind, my.ind)
                      }
                          if (length(this.ind)<1) {
                              stop("The samplePrefix value you entered doesn't",
                              " match any columns. Please check.")
                          }
                          my.counts <- my.tsrs[,this.ind]
                          my.ma <- as.matrix(my.counts)
                  }
              }

              my.names <- names(my.tsrs)
              my.cnames <- my.names[this.ind]
              n.len <- length(my.cnames)
              colData <- DataFrame(Treatment=my.cnames,
                                       row.names=my.cnames)
              my.gr <- makeGRangesFromDataFrame(my.tsrs,
                             seqnames.field="seq",                                                
                             keep.extra.columns=FALSE,
                             ignore.strand=FALSE,
                             start.field="start",
                             end.field="end",
                             strand.field="strand")
              my.se <- SummarizedExperiment(assays=list(counts=my.ma),
                                                rowRanges=my.gr,
                                            colData=colData)
              return(my.se)    
              
              message("------------------------------------------------------\n")
              message(" Done.\n")
          }
          )
