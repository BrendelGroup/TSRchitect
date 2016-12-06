#' Sets names and replicate information for experimental samples in a tssExp object
#' @param expName an S4 object of class tssExp that contains information about the experiment
#' @param sample.names unique labels (of class character) for each TSS sample within the experiment
#' @param replicate.IDs identifiers indicating which samples are biological replicates
#' @return names and replicate information for experimental samples assigned to your tssExp object
#' @export

setGeneric(
    name="setSampleID",
    def=function(expName, sample.names, replicate.IDs) {
        standardGeneric("setSampleID")
    }
    )

setMethod("setSampleID",
          signature(expName="tssExp", sample.names="character", replicate.IDs="numeric"),
          function(expName, sample.names, replicate.IDs) {
              object.name <- deparse(substitute(expName))
              exp.len <- length(expName@fileNames)

              message("... setSampleID ...")
              if (exp.len!=length(sample.names) || exp.len!=length(replicate.IDs)) {
                  stop("\nThe number of sample names and replicate IDs must be equal to the number of file names in your tssExp object.")
              }

              if (length(sample.names)!=(length(replicate.IDs))) {
                  stop("\nsample.names and replicate.IDs must have equal lengths.")
              }
              unique(sample.names) -> s.uni

              if (length(s.uni)<length(sample.names)) {
                  stop("\nEach sample name must have a unique name.")
              }

              expName@sampleNames <- sample.names
              expName@replicateIDs <- replicate.IDs
              exp.len <- length(replicate.IDs)
              rep.list <- vector(mode="list", length=exp.len)
              expName@tsrData <- rep.list

              cat("\nNames and replicate IDs were successfully assigned to tssExp\nobject \"", object.name, "\".\n\n")

              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, parent.frame())
              message(" Done.\n")
          }
          )
