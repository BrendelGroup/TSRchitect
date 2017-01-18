#' Sets names and replicate information for experimental samples in a tssObject object
#' @param experimentName an S4 object of class tssObject that contains information about the experiment.
#' @param sample.names unique labels (of class character) for each TSS sample within the experiment. Please note that importBam attaches bamData in ascending alphanumeric order, so sample.names must be arranged in this order also so that they directly correspond to the intended file.
#' @param replicate.IDs identifiers indicating which samples are biological replicates. Please note that importBam ataches bamData in ascending alphanumeric order, so replicate.IDs must be arranged in this order also so that they directly correspond to the intended file.
#' @return names and replicate information for experimental samples assigned to your tssObject object
#' @export

setGeneric(
    name="setSampleID",
    def=function(experimentName, sample.names, replicate.IDs) {
        standardGeneric("setSampleID")
    }
    )

setMethod("setSampleID",
          signature(experimentName="tssObject", sample.names="character", replicate.IDs="numeric"),
          function(experimentName, sample.names, replicate.IDs) {
              object.name <- deparse(substitute(experimentName))
              exp.len <- length(experimentName@fileNames)

              message("... setSampleID ...")
              if (exp.len!=length(sample.names) || exp.len!=length(replicate.IDs)) {
                  stop("\nThe number of sample names and replicate IDs must be equal to the number of file names in your tssObject object.")
              }

              if (length(sample.names)!=(length(replicate.IDs))) {
                  stop("\nsample.names and replicate.IDs must have equal lengths.")
              }
              unique(sample.names) -> s.uni

              if (length(s.uni)<length(sample.names)) {
                  stop("\nEach sample name must have a unique name.")
              }

              experimentName@sampleNames <- sample.names
              experimentName@replicateIDs <- replicate.IDs
              exp.len <- length(replicate.IDs)
              rep.list <- vector(mode="list", length=exp.len)
              experimentName@tsrData <- rep.list

              cat("\nNames and replicate IDs were successfully assigned to tssObject\nobject \"", object.name, "\".\n\n")

              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
