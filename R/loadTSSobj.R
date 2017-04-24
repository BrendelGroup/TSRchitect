#' @title \strong{loadTSSobj}
#' @description \code{loadTSSobj} processes alignment files in .bam format
#' from the local directory supplied.
#'
#' @param experimentTitle a descriptive title for the experiment (character).
#' @param inputDir path to the directory containing the alignment files (in
#' either .bam or .bed format, depending on your selection) (character).
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bam in \emph{inputDir} will be imported with this function.
#' @param inputType specifies the format of the TSS dataset to be imported.
#' Either "bam" or "bed" are currently accepted. (character)
#' @param isPairedEnd specifies whether the TSS profiling experiment is
#' paired-end (if TRUE) or single-end (if FALSE) (logical)
#' @param sampleNames unique labels of class character for each TSS sample
#' within the experiment (character).
#' @param replicateIDs identifiers indicating which samples are biological
#' replicates. Note that \code{loadTSSobj} imports
#' alignment data in ascending alphanumeric order, so the arguments to
#' replicateIDs must be arranged in this order also so that they directly
#' correspond to the intended file (numeric).
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
#' test.Obj <- loadTSSobj(experimentTitle="Code example", inputDir=extdata.dir,
#' inputType="bam", isPairedEnd=TRUE, sampleNames=c("sample1-rep1",
#' "sample1-rep2", "sample2-rep1","sample2-rep2"), replicateIDs=c(1,1,2,2))
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' @note All .bam files found in \emph{inputDir} will be
#' retrieved and written in ascending alphanumeric order to the
#' \emph{@fileNames} slot on the \emph{tssObject} that is created.
#'
#' @export
#' @rdname loadTSSobj-methods


setGeneric("loadTSSobj",
    function(experimentTitle, inputDir, inputType, isPairedEnd, sampleNames,
             replicateIDs)
    standardGeneric("loadTSSobj")
)

#' @rdname loadTSSobj-methods

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
                  bamFlags <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
                              isFirstMateRead=TRUE, hasUnmappedMate=FALSE,
                              isUnmappedQuery=FALSE,
                              isSecondaryAlignment=FALSE)
                  message("\nTSS data were specified to be paired-end",
                      " read alignments.")
                  myFields <- c("rname","flag","pos","qwidth","mapq",
                    "cigar","isize")
              }
              else {
                  message("\nImporting single-end reads ...\n")
                  bamFlags <- scanBamFlag(isPaired=FALSE,
                                         isUnmappedQuery=FALSE,
                                         isSecondaryAlignment=FALSE)
                  message("\nTSS data were specified to be",
                      " single-end read alignments.\n")
                  myFields <- c("rname","flag","pos",
                    "qwidth","mapq","cigar")
              }

              my.param <- ScanBamParam(flag=bamFlags, what=myFields)
              bam.paths <- tssObj@fileNames
              bv_obj <- BamViews(bam.paths)
              bv_files <- dimnames(bv_obj)[[2]]
              n.bams <- length(bv_files)
              message("\nBeginning import of ", n.bams, " bam files ...\n")
              bams.GA <- bplapply(bam.paths, readGAlignments,
                         BPPARAM = MulticoreParam(), param=my.param)
              tssObj@bamData <- bams.GA
              message("Done. Alignment data from ", n.bams,
                  " bam files have been attached to the tssObject.\n")
              message("-----------------------------------------------------\n")
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
              s.uni <- unique(sampleNames)
              if (length(s.uni)<length(sampleNames)) {
                  stop("\nEach sample name must be unique.")
              }

              tssObj@sampleNames <- sampleNames
              tssObj@replicateIDs <- replicateIDs
              exp.len <- length(replicateIDs)
              rep.list <- vector(mode="list", length=exp.len)
              tssObj@tsrData <- rep.list

              message("\nNames and replicate IDs were successfully added",
                  " to the tssObject.\n\n")
              message("-----------------------------------------------------\n")
              message(" Done.\n")
              return(tssObj)
          }
          )
