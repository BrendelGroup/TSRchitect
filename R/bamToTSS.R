#' @title \strong{bamToTSS}
#' @description \code{bamToTSS} extracts TSS information from each
#' attached .bam file in a tssObject object
#'
#' @param experimentName an S4 object of class tssObject with bam files loaded
#'
#' @return creates a list of TSSs in class \linkS4class{GRanges} for each
#' .bam file contained within \emph{experimentName} and places them in
#' the returned \emph{tssObject}.
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
                      "to TSS data conversion ...\n\n")
              }

              bam.len <- length(experimentName@bamData)
              bam.vec <- vector(mode="list", length=bam.len)

              for (i in 1:bam.len) {
                  message("Retrieving data from bam file #", i, "...\n\n")
                  experimentName@bamData[[i]] -> bam.data
                  as(bam.data,"data.frame") -> bam.df
                  bam.df[bam.df$strand=="+",] -> df.plus
                  bam.df[bam.df$strand=="-",] -> df.minus
                        gr1 <- GRanges(seqnames=df.plus$seqnames,
                                      ranges = IRanges(
                                          start=df.plus$start,
                                          end=df.plus$start
                                          ),
                                      strand=df.plus$strand
                                      )
                        gr2 <- GRanges(seqnames=df.minus$seqnames,
                                       ranges = IRanges(
                                           start=df.minus$end,
                                           end=df.minus$end
                                           ),
                                       strand=df.minus$strand
                                       )
                        c(gr1,gr2) -> gr.combined
                        sortSeqlevels(gr.combined) -> gr.combined
                        sort(gr.combined) -> gr.combined
                        gr.combined -> bam.vec[[i]]
              }

              GR.list <- GRangesList(bam.vec)
              experimentName@tssTagData <- GR.list
              experimentName@tssCountData <- vector(mode="list", length=bam.len)
              message("Done. TSS data from ", bam.len, " separate bam files" ,
                  "have been successfully\nadded to the tssObject.\n\n")
              message("---------------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
