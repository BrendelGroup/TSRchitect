#' @title \strong{initializeExp}
#' @description \code{initializeExp} creates an initialized tssObject
#' in the user's workspace
#' 
#' @param expTitle A descriptive title for the experiment (character)
#' @param experimentName The name for the \emph{tssObject} to be created
#' (character).
#' Your choice for \emph{experimentName} must lack spaces.
#' 
#' @param isPairedEnd specifies whether the TSS profiling experiment is
#' paired-end (if TRUE) or single-end (if FALSE) (logical)
#'
#' @return Creates a new \emph{tssObject} with the name \emph{experimentName}
#' that is written to the user's working environment.
#' 
#' @export

setGeneric(
    name="initializeExp",
    def=function(expTitle, experimentName, isPairedEnd) {
        standardGeneric("initializeExp")
    }
    )

setMethod("initializeExp",
          signature(expTitle="character", experimentName="character",
                    isPairedEnd="logical"),
          function(expTitle, experimentName, isPairedEnd) {
              tssObj <- new("tssObject")

              message("... initializeExp ...")

              if (isPairedEnd==TRUE) {
                  tssObj@dataType <- c("pairedEnd")
              }
              else {
                  tssObj@dataType <- c("singleEnd")
              }

              tssObj@title <- expTitle
              cat("\nThe tssObject object \"", experimentName,
                  "\" has been initialized in your workspace.\n")
              cat("-------------------------------------------------------\n")
              assign(experimentName, tssObj, parent.frame())
              message(" Done.\n")
          }
          )
