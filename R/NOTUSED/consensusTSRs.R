#' Generates a set of consensus TSRs from all conditions within the experiments loaded
#' @param experimentName - a S4 object of class tssObject that contains information about the experiment
#' @return a data frame of TSRs from all conditions
#' @export

setGeneric(
    name="consensusTSRs",
    def=function(experimentName) {
        standardGeneric("consensusTSRs")
    }
    )

setMethod("consensusTSRs",
          signature(experimentName="tssObject"),
          function(experimentName) {
              experimentName.chr <- deparse(substitute(experimentName))
              exp.type <- experimentName@dataType

              message("... consensusTSRs ...")

              ### enter rest of function here

          }
          )
