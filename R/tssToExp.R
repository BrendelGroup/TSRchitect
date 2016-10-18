#' Creates a data frame with a list of TSS positions and the corresponding counts from each
#' @param expName an S4 object of class tssExp that contains information about the experiment (including tssData)
#' @return Assigns a data frame to the expData slot of your tssExp object
#' @export 

setGeneric(
    name="tssToExp",
    def=function(expName) {
        standardGeneric("tssToExp")
    }
    )

setMethod("tssToExp",
          signature(expName="tssExp"),
          function(expName) {

              object.name <- deparse(substitute(expName))

              message("\nTSS conversion is underway.")

              if (length(expName@tssData) == 0) {

                  stop("\nAll tssData slots are empty.\nYou must load your TSS data before proceeding with this function.")
              }

              n.slots <- length(expName@tssData)

              tss.list <- vector(mode="list")

              for (i in 1:n.slots) {

                  paste("TSSset", i, sep="_") -> this.name

                  acquireTSS(expName, i) -> tss.set

                  expressionCTSS(tss.set, writeDF=FALSE) -> tss.mat

                  tss.mat -> tss.list[[this.name]]

              }

              cat("\nTSS abundance data was successfully added to your tssExp object.\n")

              return(tss.list)

#             assign(object.name, expName, envir = parent.frame())              

          }
          )
