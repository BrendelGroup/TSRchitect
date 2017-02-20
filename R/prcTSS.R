#' @title \emph{prcTSS()}
#'
#' @description An internal function that Creates an expression matrix for
#' all TSSs within a given TSS experiment (in slot tssTagData)
#'
#' @param experimentName - a S4 object of class tssObject containing
#' information in slot tssTagData
#' @param tssSet - number of the dataset to be analyzed
#' @param writeTable - if set to TRUE, writes a data frame containing
#' the TSSs positions and their abundance to your workspace
#'
#' @keywords internal
#'
#' @return \emph{prcTSS} fills the slot experimentName@tssCountData[[tssSet]]
#' in the tssObject \emph{experimentName} with a matrix of unique TSS positions
#' (rows) and observed tag counts in each position for data set tssSet;
#' precisely, each entry of the matrix consists of "seq" (chr), "TSS" (num),
#' "nTSSs" (num), and "strand" (+ or -),
#' corresponding to the sequence identifier, position, tag count, and strand,
#' respectively.
#'
#' @importFrom gtools mixedsort


setGeneric("prcTSS",
    function(experimentName, tssSet, writeTable)
    standardGeneric("prcTSS")
)

setMethod("prcTSS",
          signature(experimentName="tssObject", "numeric", "logical"),

          function(experimentName, tssSet, writeTable) {
              object.name <- deparse(substitute(experimentName))

              if (tssSet>length(experimentName@replicateIDs)) {
                  stop("The value selected exceeds the toal number of samples.")
              }

              tss <- acquireTSS(experimentName, tssSet)

              outfname <- paste("TSSset-", tssSet, sep="")
              outfname <- paste(outfname, "txt", sep=".")

              if (writeTable=="TRUE") {
                  tss.df <- tagCountTSS(tss, outfname = outfname, writeDF=TRUE)
              }
              else {
                  tss.df <- tagCountTSS(tss, outfname = outfname, writeDF=FALSE)
              }

              cat("\n... the TSS expression matrix for dataset ", tssSet,
                  " has been successfully added to\ntssObject object \"",
                  object.name, "\"\n")
              cat("---------------------------------------------------------\n")
              return(tss.df)
          }
          )
