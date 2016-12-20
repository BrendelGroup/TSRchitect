#' Creates a data frame that contains the expression of each 
#' @param expName - a S4 object of class tssExp that contains information about the experiment
#' @return alignment data (in BAM format) from the tss profiling experiments assigned to your tssExp object
#' @export 

setGeneric(
    name="tsrAbundance",
    def=function(expName) {
        standardGeneric("tsrAbundance")
    }
    )

setMethod("tsrAbundance",
          signature(expName="tssExp"),
          function(expName) {
              expName.chr <- deparse(substitute(expName))
              exp.type <- expName@dataType

              message("... tsrAbundance ...")

#### fill in rest here once method for selecting consensusTSRs is complete
              
          }
          )
