#' Summarizes the data for a TSS profiling experiment and loads them to an object of class tssObjet
#' @param experimentName - a S4 object of class tssObject that contains information about the experiment
#' @return a tabular summary of the tss profiling experiments assigned to your tssObject
#' @export

setGeneric(
    name="summarizeLibrary",
    def=function(experimentName) {
        standardGeneric("summarizeLibrary")
    }
    )

setMethod("summarizeLibrary",
          signature(experimentName="tssObject"),
          function(experimentName) {
### fill in info here
          }
          )
