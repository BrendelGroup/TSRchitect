#' @title \strong{mergeSampleData}
#' @description \code{mergeSampleData} combines samples from multiple TSS experiments into a single \linkS4class{GRanges} object
#' @param experimentName an S4 object of class \emph{tssObject} that contains information about the experiment.
#' @return tssCountData datasets will be merged (according to the \emph{sampleIDs}) and assigned to your \emph{tssObject}.
#' @importFrom GenomicRanges as.data.frame
#' @importFrom gtools mixedorder
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' mergeSampleData(experimentName=tssObjectExample)
#' @note An example similar to the one provided can be found in \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd).
#' @export

setGeneric(
    name="mergeSampleData",
    def=function(experimentName) {
        standardGeneric("mergeSampleData")
    }
    )

setMethod("mergeSampleData",
          signature(experimentName="tssObject"),
          function(experimentName) {
              object.name <- deparse(substitute(experimentName))

              message("... mergeSampleData ...")
              if (length(experimentName@tssCountData)==0) {
                  stop("\nThe slot @tssCountData is empty. Please run processTSS before proceeding with this command.\n")
              }

              if (length(experimentName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssObject object is empty. Please add sampleNames to the object.\n")
              }

              if (length(experimentName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssObject object is empty. Please add replicateIDs to the object.\n")
              }

              rep.ids <- experimentName@replicateIDs
              uni.ids <- unique(rep.ids)
              uni.ids <- uni.ids[uni.ids > 0]	# ... ignore samples with replicateID equal to zero
              exp.data <- experimentName@tssCountData
              exp.list <- vector(mode="list")

              for (i in seq_along(uni.ids)) {
                  data.frame() -> my.df
                  which(rep.ids==i) -> my.ind
                  exp.data[my.ind] -> replicate.set
                  for (j in 1:length(replicate.set)) {
                      rbind(my.df, replicate.set[[j]]) -> my.df
                  }
                  my.df <- my.df[with(my.df, order(seq, TSS)),]
                  my.df <- my.df[with(my.df, mixedorder(seq)),]
                  my.df -> exp.list[[i]]
              }

#VB: The following few lines merge the merged tssCountData into the last experimentName@tssCountDataMerged slot,
#    representing the entire collection of TSS tag counts in the experiment
 
              data.frame() -> my.df
              for (i in seq_along(uni.ids)) {
                  rbind(my.df, exp.list[[i]]) -> my.df
              }
              my.df <- my.df[with(my.df, order(seq, TSS)),]
              my.df <- my.df[with(my.df, mixedorder(seq)),]
              my.df -> exp.list[[i+1]]

              experimentName@tssCountDataMerged <- exp.list
              cat("\n... the TSS expression data has been successfully merged and added to\ntssObject object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, envir = parent.frame())
              message(" Done.\n")
          }
          )
