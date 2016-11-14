#' tssExpr
#' Creates an expression matrix for all TSSs within a given TSS experiment (in tssData)
#' @param expName - a S4 object of class tssExp containing information in slot tssData
#' @param tssNum - number of the dataset to be analyzed
#' @param writeTable if TRUE, writes a data frame containing the TSSs positions and their abundance to your workspace
#' @importFrom gtools mixedsort
#' @return creates a data frame containing tss expression for each CTSS in slot 'tssExpr' on your tssExp object
#' @export 

setGeneric(
           name="tssExpr",
           def=function(expName, tssNum, writeTable) {
               standardGeneric("tssExpr")
    }
    )

setMethod("tssExpr",
          signature(expName="tssExp", "numeric", "logical"),

          function(expName, tssNum, writeTable) {
              object.name <- deparse(substitute(expName))

              message("... tssExpr ...")
              if (tssNum>length(expName@tssData)) {
                  stop("The value selected exceeds the number of slots in tssData.")
              }
              
              tss <- acquireTSS(expName, tssNum)

              df.name <- paste("CTSSset-", tssNum, sep="")
              df.name <- paste(df.name, "txt", sep=".")              

              if (writeTable=="TRUE") {
                  tss.mat <- expressionCTSS(tss, dfName=df.name, writeDF=TRUE)
              }
              else {
                  tss.mat <- expressionCTSS(tss, dfName="my.df", writeDF=FALSE)
              }

              expName@expData[[tssNum]] <- tss.mat
              cat("\n... the TSS expression matrix for dataset ", tssNum, " has been successfully added to\ntssExp object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, envir = parent.frame())              
              message(" Done.\n")
          }
          )
