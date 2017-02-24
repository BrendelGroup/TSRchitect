#' @title \strong{loadTSSobj}
#' @description \code{loadTSSobj} processes alignment files in .bam format
#' from the local directory supplied.
#'
#' @param experimentName an S4 object of class \emph{tssObject} that
#' contains information about the experiment
#' @param inputDir path to the directory containing the alignment files in .bam
#' format (character).
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bam in \emph{inputDir} will be imported with this function.
#'
#' @return \emph{loadTSSobj} fills the slot \emph{bamData} on the returned
#' \emph{tssObject} with \linkS4class{GAlignments} objects from the
#' \bold{GenomicAlignments} package, one for each attached .bam file
#' on the \emph{fileNames} slot.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' extdata.dir <- system.file("extdata", package="TSRchitect")
#' loadTSSobj(experimentName=tssObjectExample, inputDir=extdata.dir)
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' @note All .bam files found in \emph{inputDir} will be
#' retrieved and written in ascending alphanumeric order to the
#' \emph{@fileNames} slot on the \emph{tssObject} that is created.
#'
#' @export


setGeneric("loadTSSobj",
    function(experimentTitle, inputDir, inputType, isPairedEnd, sampleNames,
             replicateIDs)
    standardGeneric("loadTSSobj")
)

setMethod("loadTSSobj",
          signature(experimentTitle="character", inputDir="character",
                    inputType="character", isPairedEnd="logical",
                    sampleNames="character", replicateIDs="numeric"),
          function(experimentTitle, inputDir, inputType, isPairedEnd,
                   sampleNames, replicateIDs) {
              
              message("... loadTSSobj ...")
              tssObj <- new("tssObject")

              tssObj@title <- experimentTitle
              tssObj.inputDir <- inputDir
              tssObj.inputType <- inputType
              if (isPairedEnd==TRUE) {
                  tssObj@dataType <- c("pairedEnd")
              }
              else {
                  tssObj@dataType <- c("singleEnd")
              }

if (inputType=="bam") {
              tss_files <- list.files(inputDir, pattern="\\.bam$",
                                      all.files=FALSE, full.names=TRUE)
              if (length(tss_files) < 1) {
                  stop("There are no .bam files in the directory you",
                       "specified, or the directory itself does not exist.",
                       "\n Please check your input for the argument 'inputDir'.")
              }
              tssObj@fileNames <- tss_files

              if(tssObj@dataType=="pairedEnd") {
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
              bam.paths <- tssObj@fileNames
              bv_obj <- BamViews(bam.paths)
              bv_files <- dimnames(bv_obj)[[2]]
              n.bams <- length(bv_files)
              cat("\nBeginning import of ", n.bams, " bam files ...\n")
              bams.GA <- bplapply(bam.paths, readGAlignments,
                         BPPARAM = MulticoreParam(), param=my.param)
              tssObj@bamData <- bams.GA
              cat("Done. Alignment data from ", n.bams,
<<<<<<< HEAD
                  " bam files have been attached to the tssObject.\n")
=======
                  " bam files have been attached to tssObject.\n")
>>>>>>> 7b41aea37ec0bd1844c567ba38f7db1d2878f3ee
              cat("---------------------------------------------------------\n")
}
if (inputType=="bed") {
stop("\nNot yet supported.  Visit again soon.\n\n")
}

              if (length(sampleNames)!=length(tssObj@fileNames)) {
                  stop("\nNumber of sampleNames must be equal to",
                       " number of input files.")
              }
              if (length(sampleNames)!=(length(replicateIDs))) {
                  stop("\nsampleNames and replicateIDs must have",
                       " equal lengths.")
              }
              unique(sampleNames) -> s.uni
              if (length(s.uni)<length(sampleNames)) {
                  stop("\nEach sample name must have a unique name.")
              }

              tssObj@sampleNames <- sampleNames
              tssObj@replicateIDs <- replicateIDs
              exp.len <- length(replicateIDs)
              rep.list <- vector(mode="list", length=exp.len)
              tssObj@tsrData <- rep.list

              cat("\nNames and replicate IDs were successfully added",
<<<<<<< HEAD
                  "to the tssObject.\n\n")
=======
                  "to tssObject.\n\n")
>>>>>>> 7b41aea37ec0bd1844c567ba38f7db1d2878f3ee
              cat("---------------------------------------------------------\n")
              message(" Done.\n")
              return(tssObj)
          }
          )
