#' @title \strong{bamToTSS}
#' @description \code{bamToTSS} extracts TSS information from each
#' attached .bam file in a tssObject object
#'
#' @param experimentName an S4 object of class tssObject with bam files loaded
#'
#' @return produces a \linkS4class{GRangesList} containing separate
#' \linkS4class{GRanges} objects for each .bam file contained within
#' \emph{experimentName}, placing them them in the returned \emph{tssObject}.
#'
#' @import BiocGenerics
#' @import methods
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' tssObjectExample <- bamToTSS(experimentName=tssObjectExample)
#'
#' @note An example similar to the one provided can be
#' found in the vignette (/inst/doc/TSRchitect.Rmd).
#'
#' @export


setGeneric("bamToTSS",
    function(experimentName)
    standardGeneric("bamToTSS")
)

setMethod("bamToTSS",
          signature(experimentName="tssObject"),
          function(experimentName) {

              message("... bamToTSS ...")
              if (length(experimentName@bamData) == 0) {
                  stop("@bamData is empty.\n\n",
                       "Please load alignment files to your tssObject.")
              }
              else {
                  message("\nBeginning .bam read alignment",
                      " to TSS data conversion ...\n\n")
              }

              bam.len <- length(experimentName@bamData)
              bam.vec <- vector(mode="list", length=bam.len)
              
              bam.df <- lapply(experimentName@bamData, as.data.frame)
              bam.gr <- lapply(bam.df, makeGRangesFromDataFrame, keep.extra.columns=FALSE)
              
              for (i in 1:bam.len) {
                  message("Retrieving data from bam file #", i, "...\n\n")
                  this.gr <- bam.gr[[i]]
                  gr.list <- S4Vectors::split(this.gr, strand(this.gr))
                  gr.plus <- gr.list$'+'
                  gr.minus <- gr.list$'-'
                  gr1 <- GRanges(seqnames=seqnames(gr.plus),
                                 ranges=IRanges(
                                     start=start(gr.plus),
                                     end=start(gr.plus)
                                     ),
                                 strand=strand(gr.plus)
                                 )
                  gr2 <- GRanges(seqnames=seqnames(gr.minus),
                                 ranges=IRanges(
                                     start=end(gr.minus),
                                     end=end(gr.minus)
                                     ),
                                 strand=strand(gr.minus)
                                 )
                        gr.combined <- c(gr1,gr2)
                        gr.combined <- sortSeqlevels(gr.combined)
                        gr.combined <- sort(gr.combined)
                        bam.vec[[i]] <- gr.combined
              }

              GR.list <- GRangesList(bam.vec)
              experimentName@tssTagData <- GR.list
              experimentName@tssCountData <- vector(mode="list", length=bam.len)
              message("Done. TSS data from ", bam.len, " separate bam files" ,
                  " have been successfully\nadded to the tssObject.\n\n")
              message("----------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
