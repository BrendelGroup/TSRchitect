#' Initializes the TSS profiling experiment
#' @param expTitle    - title for the experiment (character)
#' @param expName     - name for the tssExp object to be created in your working environment (character)
#' @param expDir      - path to the directory with the alignment files in bam format (character)
#' @param isPairedEnd - TRUE/FALSE according to whether the TSS profiling experiment was paired-end or not (logical)
#' @return Creates a tssExp object in the user's current workspace.
#' @export

setGeneric(
    name="initializeExp",
    def=function(expTitle, expName, expDir, isPairedEnd) {
        standardGeneric("initializeExp")
    }
    )

setMethod("initializeExp",
          signature(expTitle="character", expName="character", expDir="character", isPairedEnd="logical"),
          function(expTitle, expName, expDir, isPairedEnd) {
              tssObj <- new("tssExp")

              message("... initializeExp ...")
              if (expTitle == "") {
                  stop("Argument 'expTitle' of initializeExp() is empty.\n  Please provide a title for your experiment.")
              }

              if (expName == "") {
                  stop("Argument 'expName' of initializeExp() is empty.\n  Please specify a name for the tssExp object.")
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
              cat("\nThe tssExp object \"", expName, "\" has been initialized in your workspace.\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(expName, tssObj, parent.frame())
              message(" Done.\n")
          }
          )
