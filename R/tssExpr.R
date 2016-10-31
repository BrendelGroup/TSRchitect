#' tssExpr
#' Creates an expression matrix for all TSSs within a given TSS experiment (in tssData)
#' @param expName an object of tssExp format containing information in slot tssData
#' @param tssNum the number of the dataset to be analyzed
#' @param nTSSs the number of TSSs required at a given position
#' @param writeTable if TRUE, writes a data frame containing the TSSs positions and their abundance to your workspace
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

              if (tssNum>length(expName@tssData)) {
                  stop("The value selected exceeds the number of slots in tssData.")
              }
              
              tss <- acquireTSS(expName, tssNum)
              message("\nCreating expression matrix for dataset ", tssNum, "...\n")
              
              if (writeTable=="TRUE") {
              df.name <- paste("CTSS", tssNum, sep="")
              df.name <- paste(df.name, "txt", sep=".")
              tss.mat <- expressionCTSS(tss, dfName=df.name, writeDF=TRUE)
              }

              else {
              tss.mat <- expressionCTSS(tss, dfName="my.df", writeDF=FALSE)
              }

              expName@expData[[tssNum]] <- tss.mat
              message("\nA TSS expression matrix was successfully added to your tssExp object.\n")
              assign(object.name, expName, envir = parent.frame())              
          }
          )
