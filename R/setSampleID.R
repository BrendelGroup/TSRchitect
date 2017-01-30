#' @title \strong{bamToTSS}
#' @description \code{bamToTSS} Sets names and replicate information
#' for experimental samples in a \emph{tssObject} S4 object.
#' 
#' @param experimentName an S4 object of class \emph{tssObject} that
#' contains information about the experiment.
#' @param sample.names unique labels of class character for each TSS sample
#' within the experiment.
#' Please note that \code{importBam} attaches .bam data in ascending
#' alphanumeric order, so \emph{sample.names} must be arranged in this order
#' also so that they directly correspond to the intended file.
#' @param replicate.IDs identifiers indicating which samples are biological
#' replicates. As with \emph{sample.names}, note that \code{importBam} imports
#' bam data in ascending alphanumeric order, so replicate.IDs must be arranged
#' in this order also so that they directly correspond to the intended file.
#' 
#' @return names and replicate information for experimental samples assigned
#' to your \emph{tssObject} object.
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' setSampleID(experimentName=tssObjectExample,
#' sample.names=c("sample1-rep1","sample1-rep2","sample2-rep1","sample2-rep2"),
#' replicate.IDs=c(1,1,2,2))
#' #experiments 1 & 2 and 3 & 4 are replicates of samples 1 and 2, respectively
#' 
#' @note An example similar to the one provided can be found in
#' \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd)
#' 
#' @export

setGeneric(
    name="setSampleID",
    def=function(experimentName, sample.names, replicate.IDs) {
        standardGeneric("setSampleID")
    }
    )

setMethod("setSampleID",
          signature(experimentName="tssObject", sample.names="character",
                    replicate.IDs="numeric"),
          function(experimentName, sample.names, replicate.IDs) {
              object.name <- deparse(substitute(experimentName))
              exp.len <- length(experimentName@fileNames)

              message("... setSampleID ...")
              if (exp.len!=length(sample.names) || 
                  exp.len!=length(replicate.IDs)) {
                  stop("\nThe number of sample names and replicate IDs",
                       "must be equal to the number of file names in your",
                       "tssObject object.")
              }

              if (length(sample.names)!=(length(replicate.IDs))) {
                  stop("\nsample.names and replicate.IDs must have",
                       " equal lengths.")
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

              cat("\nNames and replicate IDs were successfully assigned",
                  "to tssObject\nobject \"", object.name, "\".\n\n")

              cat("-----------------------------------------------------\n")
              assign(object.name, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
