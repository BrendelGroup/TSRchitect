#' @title \emph{importBam()}
#'
#' @description \emph{importBam} processes BAM files as specified by \emph{initializeExp}.
#'
#' @param experimentName - an S4 object of class tssObject that contains information about the experiment
#'
#' @return _importBam_ fills the slot experimentName@bamData in the tssObject \emph{experimentName} with
#'         GAlignments objects from the \bold{GenomicAlignments} package, one for each parsed input
#'         BAM file.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#' @export

setGeneric(
    name="importBam",
    def=function(experimentName) {
        standardGeneric("importBam")
    }
    )

setMethod("importBam",
          signature(experimentName="tssObject"),
          function(experimentName) {
              experimentName.chr <- deparse(substitute(experimentName))
              exp.type <- experimentName@dataType

              message("... importBam ...")
              if(exp.type=="pairedEnd") {
                  message("\nImporting paired-end reads ...\n")
                  scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isFirstMateRead=TRUE, hasUnmappedMate=FALSE, isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE) -> bamFlags
                  cat("\nTSS data were specified to be paired-end read alignments.")
                  c("rname","flag","strand","pos","qwidth","mapq","cigar","isize") -> myFields
              }
              else {
                  message("\nImporting single-end reads ...\n")
                  scanBamFlag(isPaired=FALSE, isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE) -> bamFlags
                  cat("\nTSS data were specified to be single-end read alignments.\n")
                  c("rname","flag","strand","pos","qwidth","mapq","cigar") -> myFields
              }

              my.param <- ScanBamParam(flag=bamFlags, what=myFields)
              bam.paths <- experimentName@fileNames
              bv_obj <- BamViews(bam.paths)
              bv_files <- dimnames(bv_obj)[[2]]
              n.bams <- length(bv_files)
              cat("\nBeginning import of ", n.bams, " bam files ...\n")
              bams.GA <- bplapply(bam.paths, readGAlignments, BPPARAM = MulticoreParam(), param=my.param)
              experimentName@bamData <- bams.GA
              cat("Done. Alignment data from ", n.bams, " bam files have been attached to tssObject\nobject \"", experimentName.chr, "\".\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(experimentName.chr, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
