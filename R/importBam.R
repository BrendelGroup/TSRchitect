#' @title \strong{importBam}
#' @description \code{importBam} processes alignment files in .bam format
#' from the local directory supplied.
#'
#' @param experimentName an S4 object of class \emph{tssObject} that
#' contains information about the experiment
#' @param expDir path to the directory containing the alignment files in .bam
#' format (character).
#' Note that all the paths to all files in \emph{expDir} with the extension
#' .bam in \emph{expDir} will be imported with this function.
#' 
#' @return \emph{importBam} fills the slot \emph{bamData} on the
#' \emph{tssObject} with \linkS4class{GAlignments} objects from the
#' \bold{GenomicAlignments} package, one for each attached .bam file
#' on the \emph{fileNames} slot. 
#' 
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' extdata.dir <- system.file("extdata", package="TSRchitect"
#' full.names=TRUE)
#' importBam(experimentName=tssObjectExample, expDir=extdata.dir)
#'
#' @note An example similar to the one provided can be found in
#' \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd).
#' @note All .bam files found in \emph{expDir} will be
#' retrieved and written in ascending alphanumeric order to the
#' \emph{@fileNames} slot on the \emph{tssObject} that is created.
#' 
#' @export

setGeneric(
    name="importBam",
    def=function(experimentName, expDir) {
        standardGeneric("importBam")
    }
    )

setMethod("importBam",
          signature(experimentName="tssObject", expDir = "character"),
          function(experimentName, expDir) {
              object.name <- deparse(substitute(experimentName))
              exp.type <- experimentName@dataType

              message("... importBam ...")
              tss_files <- list.files(expDir, pattern="\\.bam$",
                           all.files=FALSE, full.names=TRUE)

              if (length(tss_files) < 1) {
                  stop("There are no .bam files in the directory you",
                       "specified, or the directory itself does not exist.",
                       "\n Please check your input for the argument 'expDir'.")
              }

              experimentName@fileNames <- tss_files
              message("\nLinks to .bam files have been assigned...\n")
              
              if(exp.type=="pairedEnd") {
                  message("\nImporting paired-end reads ...\n")
                  scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
                              isFirstMateRead=TRUE, hasUnmappedMate=FALSE,
                              isUnmappedQuery=FALSE,
                              isSecondaryAlignment=FALSE) -> bamFlags
                  cat("\nTSS data were specified to be paired-end",
                      "read alignments.")
                  c("rname","flag","strand","pos","qwidth","mapq",
                    "cigar","isize") -> myFields
              }
              else {
                  message("\nImporting single-end reads ...\n")
                  scanBamFlag(isPaired=FALSE, isUnmappedQuery=FALSE,
                              isSecondaryAlignment=FALSE) -> bamFlags
                  cat("\nTSS data were specified to be",
                      "single-end read alignments.\n")
                  c("rname","flag","strand","pos",
                    "qwidth","mapq","cigar") -> myFields
              }

              my.param <- ScanBamParam(flag=bamFlags, what=myFields)
              bam.paths <- experimentName@fileNames
              bv_obj <- BamViews(bam.paths)
              bv_files <- dimnames(bv_obj)[[2]]
              n.bams <- length(bv_files)
              cat("\nBeginning import of ", n.bams, " bam files ...\n")
              bams.GA <- bplapply(bam.paths, readGAlignments,
                         BPPARAM = MulticoreParam(), param=my.param)
              experimentName@bamData <- bams.GA
              cat("Done. Alignment data from ", n.bams,
                  " bam files have been attached to tssObject\nobject \"",
                  object.name, "\".\n")
              cat("-------------------------------------------------------\n")
              assign(object.name, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
