#' processTSS
#' Creates an expression matrix for all TSSs within a given TSS experiment (in tssData)
#' @param experimentName - a S4 object of class tssObject containing information in slot tssData
#' @param tssSet - number of the dataset to be analyzed
#' @param writeTable if TRUE, writes a data frame containing the TSSs positions and their abundance to your workspace
#' @importFrom gtools mixedsort
#' @return creates a data frame containing tss expression for each TSS position in slot 'expData' on your tssObject object
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
              if (tssSet>length(experimentName@tssData)) {
                  stop("The value selected exceeds the number of slots in tssData.")
              }

              tss <- acquireTSS(experimentName, tssSet)

              df.name <- paste("TSSset-", tssSet, sep="")
              df.name <- paste(df.name, "txt", sep=".")

              if (writeTable=="TRUE") {
                  tss.mat <- tagCountTSS(tss, dfName=df.name, writeDF=TRUE)
              }
              else {
                  tss.mat <- tagCountTSS(tss, dfName="my.df", writeDF=FALSE)
              }

              experimentName@expData[[tssSet]] <- tss.mat
              cat("\n... the TSS expression matrix for dataset ", tssSet, " has been successfully added to\ntssObject object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, envir = parent.frame())
              message(" Done.\n")
          }
          )
