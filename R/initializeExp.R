#' Initializes the TSS profiling experiment
#' @param expTitle    - title for the experiment (character)
#' @param experimentName     - name for the tssObject object to be created in your working environment (character)
#' @param expDir      - path to the directory with the alignment files in bam format (character)
#' @param isPairedEnd - TRUE/FALSE according to whether the TSS profiling experiment was paired-end or not (logical)
#' @return Creates a tssObject object in the user's current workspace.
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
