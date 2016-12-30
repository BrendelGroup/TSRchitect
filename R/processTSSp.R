#' processTSSp
#' Creates an expression matrix for all TSSs within a given TSS experiment (in tssTagData)
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param tssSet - number of the dataset to be analyzed
#' @param writeTable if TRUE, writes a data frame containing the TSSs positions and their abundance to your workspace
#' @importFrom gtools mixedsort
#' @return creates a data frame containing tss expression for each TSS to be stored in slot 'tssExpression' on your tssObject object
#' @export

setGeneric(
           name="processTSSp",
           def=function(experimentName, tssSet, writeTable) {
               standardGeneric("processTSSp")
    }
    )

setMethod("processTSSp",
          signature(experimentName="tssObject", "numeric", "logical"),

          function(experimentName, tssSet, writeTable) {
              object.name <- deparse(substitute(experimentName))

              message("... processTSSp ...")
              if (tssSet>length(experimentName@tssTagData)) {
                  stop("The value selected exceeds the number of slots in tssTagData.")
              }

              tss <- acquireTSS(experimentName, tssSet)

              df.name <- paste("TSSset-", tssSet, sep="")
              df.name <- paste(df.name, "txt", sep=".")

              if (writeTable=="TRUE") {
                  tss.mat <- tagCountTSS(tss, dfName=df.name, writeDF=TRUE)
              }
              else {
                  tss.mat <- tagCountTSS(tss, dfName=df.name, writeDF=FALSE)
              }

              cat("\n... the TSS expression matrix for dataset ", tssSet, " has been successfully added to\ntssObject object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              message(" Done.\n")
              return(tss.mat)
          }
          )
