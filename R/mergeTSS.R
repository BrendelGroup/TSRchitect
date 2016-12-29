#' Combines samples from different tss experiments into a single GRanges object
#' @param experimentName - an S4 object of class tssObject that contains information about the experiment.
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom S4Vectors Rle
#' @importFrom IRanges IRanges
#' @return merged tss profiling experiments according to their assigned sampleIDs to the tssObject
#' @export

setGeneric(
    name="mergeTSS",
    def=function(experimentName) {
        standardGeneric("mergeTSS")
    }
    )

setMethod("mergeTSS",
          signature(experimentName="tssObject"),
          function(experimentName) {
              experimentName.chr <- deparse(substitute(experimentName))

              if (length(experimentName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssObject is empty. Please add sampleNames to the object.\n")
              }

              if (length(experimentName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssObject is empty. Please add replicateIDs to the object.\n")
              }

              if (length(experimentName@tssData) < 2) {
                  stop("\nThere are less than two tssData slots loaded on your tssObject object, so merging is not possible.")
              }

              rep.ids <- experimentName@replicateIDs
              uni.ids <- unique(rep.ids)
              tss.data <- experimentName@tssData
              gr.list <- GRangesList()

              for (i in seq_along(uni.ids)) {
                  i -> sample.num
                  which(rep.ids==sample.num) -> my.ind
                  tss.data[my.ind] -> replicate.set
                  GRanges(seqnames=Rle(), ranges=NULL, strand=NULL, seqLengths=NULL, seqinfo=NULL) -> gr.set
                  #GRanges(seqnames=NULL, ranges=NULL, strand=NULL, seqLengths=NULL, seqinfo=NULL) -> gr.set
                  for (j in 1:length(replicate.set)) {
                      replicate.set[[j]] -> this.gr
                      c(gr.set, this.gr) -> gr.set
                      }
                  gr.set -> gr.list[[i]]
              }

              experimentName@tssDataMerged <- gr.list
              message("\nTSS data has been merged and assigned to your tssObject.\n")
              assign(experimentName.chr, experimentName, envir = parent.frame())
          }
          )
