#' @title \strong{inputToTSS}
#' @description \code{inputToTSS} extracts TSS information from each
#' attached .bam or .bed file in a tssObject object
#'
#' @param experimentName an S4 object of class tssObject with bam files loaded
#'
#' @return produces a \linkS4class{GRangesList} containing separate
#' \linkS4class{GRanges} objects for each .bam file contained within
#' \emph{experimentName}, placing them them in the returned \emph{tssObject}.
#' 
#' @import BiocGenerics
#' @importFrom GenomicRanges granges GRanges GRangesList
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors first second split
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' tssObjectExample <- inputToTSS(experimentName=tssObjectExample)
#'
#' @note An example similar to the one provided can be
#' found in the vignette (/inst/doc/TSRchitect.Rmd).
#'
#' @export
#' @rdname inputToTSS-methods


setGeneric("inputToTSS",
    function(experimentName)
    standardGeneric("inputToTSS")
)

#' @rdname inputToTSS-methods

setMethod("inputToTSS",
          signature(experimentName="tssObject"),
          function(experimentName) {
              message("... inputToTSS ...")
              if (length(experimentName@bamData)>0) {
                  message("\nBeginning input to",
                      " TSS data conversion ...\n\n")
                  bam.len <- length(experimentName@bamData)
                  bam.vec <- vector(mode="list", length=bam.len)
                  bam.df <- lapply(experimentName@bamData, as.data.frame)
                  bam.gr <- lapply(bam.df, makeGRangesFromDataFrame, keep.extra.columns=FALSE)
              for (i in 1:bam.len) {
                  message("Retrieving data from bam file #", i, "...\n\n")
                  this.gr <- bam.gr[[i]]
                  gr.list <- S4Vectors::split(this.gr, strand(this.gr))
                  gr.plus <- gr.list$'+'
                  gr.minus <- gr.list$'-'
                  gr1 <- GRanges(seqnames=seqnames(gr.plus),
                                 ranges=IRanges(
                                     start=start(gr.plus),
                                     end=start(gr.plus)
                                     ),
                                 strand=strand(gr.plus)
                                 )
                  gr2 <- GRanges(seqnames=seqnames(gr.minus),
                                 ranges=IRanges(
                                     start=end(gr.minus),
                                     end=end(gr.minus)
                                     ),
                                 strand=strand(gr.minus)
                                 )
              gr.combined <- c(gr1,gr2)
              gr.combined <- sortSeqlevels(gr.combined)
              gr.combined <- sort(gr.combined)
              bam.vec[[i]] <- gr.combined
              }
              GR.list <- GRangesList(bam.vec)
              experimentName@tssTagData <- GR.list
              experimentName@tssCountData <- vector(mode="list", length=bam.len)
              message("Done. TSS data from ", bam.len, " separate bam files" ,
                  " have been successfully\nadded to the tssObject.\n\n")
              message("----------------------------------------------------\n")
              }
           if (length(experimentName@bedData)>0) {
              message("\nBeginning input to",
                     " TSS data conversion ...\n\n")
              if (class(experimentName@bedData[[1]])=="GRanges") { #i.e bed single-end
              bed.len <- length(experimentName@bedData)
              bed.vec <- vector(mode="list", length=bed.len)
              bed.gr <- experimentName@bedData
              for (i in 1:bed.len) {
                  message("Retrieving data from bed file #", i, "...\n\n")
                  this.gr <- bed.gr[[i]]
                  gr.list <- S4Vectors::split(this.gr, strand(this.gr))
                  gr.plus <- gr.list$'+'
                  gr.minus <- gr.list$'-'
                  gr1 <- GRanges(seqnames=seqnames(gr.plus),
                                 ranges=IRanges(
                                     start=start(gr.plus),
                                     end=start(gr.plus)
                                     ),
                                 strand=strand(gr.plus)
                                 )
                  gr2 <- GRanges(seqnames=seqnames(gr.minus),
                                 ranges=IRanges(
                                     start=end(gr.minus),
                                     end=end(gr.minus)
                                     ),
                                 strand=strand(gr.minus)
                                 )
                  gr.combined <- c(gr1,gr2)
                  gr.combined <- sortSeqlevels(gr.combined)
                  gr.combined <- sort(gr.combined)
                  bed.vec[[i]] <- gr.combined
              }
                  if (exists("bam.len")==TRUE) {
                      combined.len <- bam.len+bed.len
                      start.pos <- bam.len+1
                      combined.vec <- vector(mode="list", length=combined.len)
                      combined.vec[1:bam.len] <- bam.vec
                      combined.vec[start.pos:combined.len] <- bed.vec
                      GR.list <- GRangesList(combined.vec)
                      experimentName@tssTagData <- GR.list
                      experimentName@tssCountData <- vector(mode="list", length=combined.len)
                  }
                  else {
                     GR.list <- GRangesList(bed.vec)
                     experimentName@tssTagData <- GR.list
                     experimentName@tssCountData <- vector(mode="list", length=bed.len)
                 }
              GR.list <- GRangesList(bed.vec)
              experimentName@tssTagData <- GR.list
              experimentName@tssCountData <- vector(mode="list", length=bed.len)
              message("Done. TSS data from ", bed.len, " separate bed files" ,
                  " have been successfully\nadded to the tssObject.\n\n")
              message("----------------------------------------------------\n")
          }
              if (class(experimentName@bedData[[1]])=="Pairs") { #i.e. bed paired-end (bedpe)
              bedpe.len <- length(experimentName@bedData)
              bedpe.vec <- vector(mode="list", length=bedpe.len)
              bed.gr <- lapply(experimentName@bedData, first)
              for (i in 1:bedpe.len) {
                  message("Retrieving data from bedpe file #", i, "...\n\n")
                  this.gr <- bed.gr[[i]]
                  gr.list <- S4Vectors::split(this.gr, strand(this.gr))
                  gr.plus <- gr.list$'+'
                  gr.minus <- gr.list$'-'
                  gr1 <- GRanges(seqnames=seqnames(gr.plus),
                                 ranges=IRanges(
                                     start=start(gr.plus),
                                     end=start(gr.plus)
                                     ),
                                 strand=strand(gr.plus)
                                 )
                  gr2 <- GRanges(seqnames=seqnames(gr.minus),
                                 ranges=IRanges(
                                     start=end(gr.minus),
                                     end=end(gr.minus)
                                     ),
                                 strand=strand(gr.minus)
                                 )
                  gr.combined <- c(gr1,gr2)
                  gr.combined <- sortSeqlevels(gr.combined)
                  gr.combined <- sort(gr.combined)
                  bedpe.vec[[i]] <- gr.combined
              }
              GR.list <- GRangesList(bedpe.vec)
                  if (exists("bam.len")==TRUE) {
                      combined.len <- bam.len+bedpe.len
                      start.pos <- bam.len+1
                      combined.vec <- vector(mode="list", length=combined.len)
                      combined.vec[1:bam.len] <- bam.vec
                      combined.vec[start.pos:combined.len] <- bedpe.vec
                      GR.list <- GRangesList(combined.vec)
                      experimentName@tssTagData <- GR.list
                      experimentName@tssCountData <- vector(mode="list", length=combined.len)
                  }
                 else {
                     GR.list <- GRangesList(bedpe.vec)
                     experimentName@tssTagData <- GR.list
                     experimentName@tssCountData <- vector(mode="list", length=bedpe.len)
                 }
              message("Done. TSS data from ", bedpe.len, " separate bed files" ,
                  " have been successfully\nadded to the tssObject.\n\n")
              message("----------------------------------------------------\n")
          }
          }
              message(" Done.\n")
              return(experimentName)
          }
          )
