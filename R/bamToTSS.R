#' bamToTSS
#' Extracts TSS information from each bam file in a tssExp object
#' @param expName - a S4 object of class tssExp with bam files loaded
#' @return creates a list of TSSs in class GRanges for each bam file contained within expName and places them in the tssExp object
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom BiocGenerics start end
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @export

setGeneric(
           name="bamToTSS",
           def=function(expName) {
               standardGeneric("bamToTSS")
    }
    )

setMethod("bamToTSS",
          signature(expName="tssExp"),
          function(expName) {
              object.name <- deparse(substitute(expName))

              message("... bamToTSS ...")
              if (length(expName@bamData) == 0) {
                  stop("Slot @bamData is empty.\n\n Please load alignment files to your tssExp object.")
              }
              else {
                  cat("\nBeginning .bam read alignment to TSS data conversion ...\n\n")
              }

              bam.len <- length(expName@bamData)
              bam.vec <- vector(mode="list", length=bam.len)

              for (i in 1:bam.len) {
                  cat("Retrieving data from bam file #", i, "...\n\n")
                  expName@bamData[[i]] -> bam.data
                  as(bam.data,"data.frame") -> bam.df
                  bam.df[bam.df$strand=="+",] -> df.plus
                  bam.df[bam.df$strand=="-",] -> df.minus
                        gr1 <- GRanges(seqnames=df.plus$seqnames,
                                      ranges = IRanges(
                                          start=df.plus$start,
                                          end=df.plus$start
                                          ),
                                      strand=df.plus$strand
                                      )
                        gr2 <- GRanges(seqnames=df.minus$seqnames,
                                       ranges = IRanges(
                                           start=df.minus$end,
                                           end=df.minus$end
                                           ),
                                       strand=df.minus$strand
                                       )
                        c(gr1,gr2) -> gr.combined
                        sortSeqlevels(gr.combined) -> gr.combined
                        sort(gr.combined) -> gr.combined
                        gr.combined -> bam.vec[[i]]
              }

              GR.list <- GRangesList(bam.vec)
              expName@tssData <- GR.list
              expName@expData <- vector(mode="list", length=bam.len)
              cat("Done. TSS data from ", bam.len, " separate bam files have been successfully added to\ntssExp object \"", object.name, "\".\n\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, expName, envir = parent.frame())
              message(" Done.\n")
          }
          )
