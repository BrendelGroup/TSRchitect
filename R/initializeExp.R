#' @title \strong{initializeExp}
#' @description \code{initializeExp} initializes the TSS profiling experiment
#' @param expTitle A descriptive title for the experiment (character)
#' @param experimentName name for the \emph{tssObject} to be created in your working environment (character). Your choice for \emph{experimentName} must lack spaces.
#' @param expDir path to the directory with the alignment files in .bam format (character). Note that all the paths to all files in \emph{expDir} with the extension .bam in \emph{expDir} will be imported with this function.
#' @param isPairedEnd specifies whether the TSS profiling experiment is paired-end (if TRUE) or single-end (if FALSE) (logical)
#' @return Creates a new \emph{tssObject} with the name \emph{experimentName} that is written to the user's working environment.
#' @note \code{initializeExp("TSS Object Example", experimentName=tssObjectExample, expDir="extdata", isPairedEnd=TRUE)}
#' @note Please note that all .bam files found in \emph{expDir} will be retrieved and written in ascending alphanumeric order to the \emph{@fileNames} slot on the \emph{tssObject} that is created.
#'
#' @importFrom methods new
#' @export

setGeneric(
    name="initializeExp",
    def=function(expTitle, experimentName, expDir, isPairedEnd) {
        standardGeneric("initializeExp")
    }
    )

setMethod("initializeExp",
          signature(expTitle="character", experimentName="character", expDir="character", isPairedEnd="logical"),
          function(expTitle, experimentName, expDir, isPairedEnd) {
              tssObj <- new("tssObject")

              message("... initializeExp ...")
              if (expTitle == "") {
                  stop("Argument 'expTitle' of initializeExp() is empty.\n  Please provide a title for your experiment.")
              }

              if (experimentName == "") {
                  stop("Argument 'experimentName' of initializeExp() is empty.\n  Please specify a name for the tssObject object.")
              }

              if (expDir == "") {
                  stop("Argument 'expDir' of initializeExp() is empty.\n  Please supply the name of the directory containing the read alignment files (in .bam format) for this experiment.")
              }

              if (isPairedEnd==TRUE) {
                  tssObj@dataType <- c("pairedEnd")
              }
              else {
                  tssObj@dataType <- c("singleEnd")
              }

              tss_files <- list.files(expDir, pattern="\\.bam$", all.files=FALSE, full.names=TRUE)

              if (length(tss_files) < 1) {
                  stop("There are no .bam files in the directory you specified, or the directory itself does not exist.\n  Please check your input for argument 'expDir' of initializeExp().")
              }

              tssObj@title <- expTitle
              tssObj@fileNames <- tss_files
              cat("\nThe tssObject object \"", experimentName, "\" has been initialized in your workspace.\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(experimentName, tssObj, parent.frame())
              message(" Done.\n")
          }
          )
