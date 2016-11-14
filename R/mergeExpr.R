#' Combines samples from two different tss experiments into a single GRanges object
#' @param expName an S4 object of class tssExp that contains information about the experiment.
#' @return expData datasets will be merged (according to the sampleIDs) and assigned to your tssExp object
#' @importFrom GenomicRanges as.data.frame
#' @export 

setGeneric(
    name="mergeExpr",
    def=function(expName) {
        standardGeneric("mergeExpr")
    }
    )

setMethod("mergeExpr",
          signature(expName="tssExp"),
          function(expName) {
              object.name <- deparse(substitute(expName))

              if (length(expName@expData)==0) {
                  stop("\nThe slot @expData is empty. Please run tssExpr before proceeding with this command.\n")
              }

              if (length(expName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssExp object is empty. Please add sampleNames to the object.\n")
              }

              if (length(expName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssExp object is empty. Please add replicateIDs to the object.\n")
              }

              rep.ids <- expName@replicateIDs
              uni.ids <- unique(rep.ids)
              exp.data <- expName@expData
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

              expName@expDataMerged <- exp.list
              cat("\n... the TSS expression data has been successfully merged and added to\ntssExp object \"", object.name, "\"\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, envir = parent.frame())
              message(" Done.\n")
          }
          )
