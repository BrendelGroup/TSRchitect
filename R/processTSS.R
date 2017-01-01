#' @title \emph{processTSS()}
#'
#' @description \emph{processTSS} Creates an expression matrix for all TSSs within a given TSS experiment (in tssTagData)
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param tssSet - number of the dataset to be analyzed
#' @param writeTable - if set to TRUE, writes a data frame containing the TSSs positions and their abundance to your workspace
#'
#' @return \emph{processTSS} fills the slot experimentName@tssCountData[[tssSet]] in
#'         the tssObject \emph{experimentName} with a matrix of unique TSS positions
#'         (rows) and observed tag counts in each position for data set tssSet;
#'         precisely, each entry of the matrix consists of "seq" (chr), "TSS" (num), "nTSSs" (num), and "strand" (+ or -),
#          corresponding to the sequence identifier, position, tag count, and strand, respectively.
#'
#' @importFrom gtools mixedsort
#' @export

setGeneric(
           name="processTSS",
           def=function(experimentName, tssSet, writeTable) {
               standardGeneric("processTSS")
    }
    )

setMethod("processTSS",
          signature(experimentName="tssObject", "numeric", "logical"),

          function(experimentName, tssSet, writeTable) {
              object.name <- deparse(substitute(experimentName))

              message("... processTSS ...")
              if (tssSet>length(experimentName@tssTagData)) {
                  stop("The value selected exceeds the number of slots in tssTagData.")
              }

              tss <- acquireTSS(experimentName, tssSet)

              df.name <- paste("TSSset-", tssSet, sep="")
              df.name <- paste(df.name, "txt", sep=".")

              if (writeTable=="TRUE") {
                  tss.df <- tagCountTSS(tss, dfName=df.name, writeDF=TRUE)
              }
              else {
                  tss.df <- tagCountTSS(tss, dfName=df.name, writeDF=FALSE)
              }

              cat("\n... the TSS expression matrix for dataset ", tssSet, " has been successfully added to\ntssObject object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              return(tss.df)
              message(" Done.\n")
          }
          )
