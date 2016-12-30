#' Combines samples from two different tss experiments into a single GRanges object
#' @param experimentName - an S4 object of class tssObject that contains information about the experiment.
#' @return tssCountData datasets will be merged (according to the sampleIDs) and assigned to your tssObject object
#' @importFrom GenomicRanges as.data.frame
#' @export

setGeneric(
    name="mergeExpr",
    def=function(experimentName) {
        standardGeneric("mergeExpr")
    }
    )

setMethod("mergeExpr",
          signature(experimentName="tssObject"),
          function(experimentName) {
              object.name <- deparse(substitute(experimentName))

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
              exp.data <- experimentName@tssCountData
              exp.list <- vector(mode="list")

              for (i in seq_along(uni.ids)) {
                  i -> sample.num
                  which(rep.ids==sample.num) -> my.ind
                  exp.data[my.ind] -> replicate.set
                  data.frame() -> my.df
                  for (j in 1:length(replicate.set)) {
                      replicate.set[[j]] -> this.df
                      rbind(my.df, this.df) -> my.df
                      my.df <- my.df[with(my.df, order(chr, CTSS)),]
                  }
                  my.df -> exp.list[[i]]
              }

              experimentName@tssCountDataMerged <- exp.list
              cat("\n... the TSS expression data has been successfully merged and added to\ntssObject object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, envir = parent.frame())
              message(" Done.\n")
          }
          )
