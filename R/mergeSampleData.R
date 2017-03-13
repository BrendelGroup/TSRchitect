#' @title \strong{mergeSampleData}
#' @description \code{mergeSampleData} combines samples from multiple TSS
#' experiments into a single \linkS4class{GRanges} object
#'
#' @param experimentName an S4 object of class \emph{tssObject} that contains
#' information about the experiment.
#'
#' @return tssCountData datasets are merged (according to the
#' \emph{sampleIDs}) and put in the tssCountDataMerged slot in the returned
#' \emph{tssObject}.
#'
#' @importFrom GenomicRanges as.data.frame
#' @importFrom gtools mixedorder
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' tssObjectExample <- mergeSampleData(experimentName=tssObjectExample)
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' @export


setGeneric("mergeSampleData",
    function(experimentName)
    standardGeneric("mergeSampleData")
)

setMethod("mergeSampleData",
          signature(experimentName="tssObject"),
          function(experimentName) {

              message("... mergeSampleData ...")
              if (length(experimentName@tssCountData)==0) {
                  stop("\nThe slot @tssCountData is empty.",
                       "Please run processTSS before proceeding with",
                       "this command.\n")
              }

              if (length(experimentName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssObject",
                       "object is empty. Please add sampleNames to",
                       "the object.\n")
              }

              if (length(experimentName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssObject",
                       "object is empty.\n",
                       "Please add replicateIDs to the object.\n")
              }

              rep.ids <- experimentName@replicateIDs
              uni.ids <- unique(rep.ids)
              uni.ids <- uni.ids[uni.ids > 0]
# ignores samples with replicateID equal to zero
              exp.data <- experimentName@tssCountData
              exp.list <- vector(mode="list")

              for (i in seq_along(uni.ids)) {
                  my.df <- data.frame()
                  my.ind <- which(rep.ids==i)
                  replicate.set <- exp.data[my.ind]
                  for (j in 1:length(replicate.set)) {
                      my.df <- rbind(my.df, replicate.set[[j]])
                  }
                  my.df <- my.df[with(my.df, order(seq, TSS)),]
                  my.df <- my.df[with(my.df, mixedorder(seq)),]
                  exp.list[[i]] <- my.df
              }

# The following few lines merge the merged tssCountData into the last
# experimentName@tssCountDataMerged slot, representing the entire
# collection of TSS tag counts in the experiment

              my.df <- data.frame()
              for (i in seq_along(uni.ids)) {
                  my.df <- rbind(my.df, exp.list[[i]])
              }
              my.df <- my.df[with(my.df, order(seq, TSS)),]
              my.df <- my.df[with(my.df, mixedorder(seq)),]
              exp.list[[i+1]] <- my.df

              experimentName@tssCountDataMerged <- exp.list
              message("\n... the TSS expression data have been successfully merged",
                    "\nand added to the tssObject object.\n")
              message("---------------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
