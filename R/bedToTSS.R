#' @title \strong{bedToTSS}
#' @description \code{bedToTSS} extracts TSS information from each
#' attached .bed file in a tssObject object
#'
#' @param experimentName an S4 object of class tssObject with bed files loaded
#'
#' @return produces a \linkS4class{GRangesList} containing separate
#' \linkS4class{GRanges} objects for each .bed file contained within
#' \emph{experimentName}, placing them them in the returned \emph{tssObject}.
#' 
#' @import BiocGenerics
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#'
#' @note An example similar to the one provided can be
#' found in the vignette (/inst/doc/TSRchitect.Rmd).
#'
#' @export
#' @rdname bedToTSS-methods


setGeneric("bedToTSS",
    function(experimentName)
    standardGeneric("bedToTSS")
)

#' @rdname bedToTSS-methods

setMethod("bedToTSS",
          signature(experimentName="tssObject"),
          function(experimentName) {

              message("... bedToTSS ...")
              if (length(experimentName@bedData) == 0) {
                  stop("@bedData is empty.\n\n",
                       "Please load alignment files to your tssObject.")
              } else {
                  message("\nBeginning .bed file",
                      " to TSS data conversion ...\n\n")
              }

              bed.len <- length(experimentName@bedData)
              bed.gr <- experimentName@bedData
              GR.list <- GRangesList(bed.gr)
              experimentName@tssTagData <- GR.list
              experimentName@tssCountData <- vector(mode="list", length=bed.len)
              message("Done. TSS data from ", bed.len, " separate bed files" ,
                  " have been successfully\nadded to the tssObject.\n\n")
              message("----------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
