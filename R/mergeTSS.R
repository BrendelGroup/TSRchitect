#' Combines samples from two different 
#' @param expName an S4 object of class tssExp that contains information about the experiment.
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @return merged tss profiling experiments according to their assigned sampleIDs to your tssExp object
#' @export 

setGeneric(
    name="mergeTSSs",
    def=function(expName) {
        standardGeneric("mergeTSSs")
    }
    )

setMethod("mergeTSSs",
          signature(expName="tssExp"),
          function(expName) {
              expName.chr <- deparse(substitute(expName))

              if (length(expName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssExp object is empty. Please add sampleNames to the object.\n")
              }

              if (length(expName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssExp object is empty. Please add replicateIDs to the object.\n")
              }

              if (length(expName@tssData) < 2) {
                  stop("\nThere are less than two tssData slots loaded on your tssExp object, so merging is not possible.")
              }

              rep.ids <- expName@replicateIDs
              uni.ids <- unique(rep.ids)
              tss.data <- expName@tssData
              gr.list <- GRangesList()
              
              for (i in seq_along(uni.ids)) {
                  i -> sample.num
                  which(rep.ids==sample.num) -> my.ind
                  tss.data[my.ind] -> replicate.set
                  for (j in 1:length(replicate.set)) {
                      GRanges(seqnames=Rle(), ranges=NULL, strand=NULL, seqLengths=NULL, seqinfo=NULL) -> gr.set
                      replicate.set[[j]] -> this.gr
                      c(gr.set, this.gr) -> gr.set
                  }
                  gr.set -> gr.list[[i]]
              }

              expName@tssDataMerged <- gr.list 
              message("\nTSS data has been merged and assigned to your tssExp object.\n")
              assign(expName.chr, expName, envir = parent.frame())
          }
          )
