#' Initializes the TSS profiling experiment
#' @param expTitle a title for the experiment (character)
#' @param objName a name the object to be created in your working environment (character)
#' @param expDir a path to the directory with the alignment files in bam format (character)
#' @param isPairedEnd a logical indicating whether the TSS profiling experiment was paired-end (logical)
#' @return Creates a tssExp object in the user's current workspace.
#' @export 

setGeneric(
    name="initializeExp",
    def=function(expTitle, objName, expDir, isPairedEnd) {
        standardGeneric("initializeExp")
    }
    )

setMethod("initializeExp",
          signature(expTitle="character", objName="character", expDir="character", isPairedEnd="logical"),
          function(expTitle, objName, expDir, isPairedEnd) {
              tssObj <- new("tssExp")
              
              if (is.na(expTitle)) {
                  stop("Argument 'expTitle' is empty. Please supply a title for your experiment.")
              }

              if (is.na(objName)) {
                  stop("Argument 'objName' is empty. Please supply a name for the tssExp object.")
                   }

              if (is.na(expDir)) {
                stop("Argument 'expDir' is empty. Please supply the full paths location of the alignment files (BAMs) for this experiment.")
             }
              
              if (isPairedEnd==TRUE) {
                  tssObj@dataType <- c("pairedEnd")
              }

              else {
                  tssObj@dataType <- c("singleEnd")
              }

              tss_files <- list.files(expDir, pattern="\\.bam$", all.files=FALSE,full.names=TRUE)

              if (length(tss_files) < 1) {
                  stop("There are no .bam files in the directory you supplied, or the directory itself does not exist.\n Please correct your input for 'expDir'.")
              }

              tssObj@title <- expTitle
              tssObj@fileNames <- tss_files
              message("The tssExp object:", objName, "is now in your workspace. \n")
              assign(objName, tssObj, parent.frame()) 
          }
          )

              
  
