#' Generates a set of consensus TSRs from all conditions within the experiments loaded
#' @param expName - a S4 object of class tssExp that contains information about the experiment
#' @return a data frame of TSRs from all conditions
#' @export

setGeneric(
    name="consensusTSRs",
    def=function(expName) {
        standardGeneric("consensusTSRs")
    }
    )

setMethod("consensusTSRs",
          signature(expName="tssExp"),
          function(expName) {
              expName.chr <- deparse(substitute(expName))
              exp.type <- expName@dataType

              message("... consensusTSRs ...")

              ### enter rest of function here

          }
          )
