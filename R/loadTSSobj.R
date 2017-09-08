#' @title \strong{loadTSSobj}
#' @description \code{loadTSSobj} processes alignment files in .bam format
#' from the local directory supplied.
#'
#' @param experimentTitle a descriptive title for the experiment (character).
#' @param inputDirBAM path to the directory containing the alignment files (in
#' either .bam format) (character). 
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bam in \emph{inputDirBAM} will be imported with this function.
#' @param inputDirBED path to the directory containing the input files (in
#' BED format) (character).
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bed in \emph{inputDirBED} will be imported with this function.
#' @param isPairedBAM if the input is in BAM format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE) (logical)
#' @param isPairedBED if the input is in BED format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE) (logical). Note: if TRUE, the input data must be in
#' bedpe format, as described here:
#' http://bedtools.readthedocs.io/en/latest/content/general-usage.html
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
#' @importFrom rtracklayer import.bed
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
           function(experimentTitle, inputDirBAM, inputDirBED,
                    isPairedBAM, isPairedBED, sampleNames,
                    replicateIDs)
    standardGeneric("loadTSSobj")
)

#' @rdname loadTSSobj-methods

setMethod("loadTSSobj", #both BAM and BED files
          signature(experimentTitle="character", inputDirBAM="ANY",
                    inputDirBED="ANY", isPairedBAM="ANY",
                    isPairedBED="ANY", sampleNames="character",
                    replicateIDs="numeric"),
          function(experimentTitle, inputDirBAM, inputDirBED=NA,
                   isPairedBAM, isPairedBED=NA, sampleNames,
                   replicateIDs) {
              
              message("... loadTSSobj ...")
              tssObj <- new("tssObject")

              tssObj@title <- experimentTitle

              if (is.character(inputDirBAM)) {
              tssObj@inputDirBAM <- inputDirBAM
             }
              if (is.character(inputDirBED)) {
              tssObj@inputDirBED <- inputDirBED
            }
              if ((is.character(inputDirBAM)==FALSE) & (is.character(inputDirBED)==FALSE)) {
                  stop("No input valide directories have been supplied.")
              }
              if ((is.character(inputDirBAM)==TRUE) & (is.logical(isPairedBAM)==FALSE)) {
                  stop("Arguments for inputDirBAM and isPairedBAM must be supplied.")
              }
              if ((is.character(inputDirBED)==TRUE) & (is.logical(isPairedBED)==FALSE)) {
                  stop("Arguments for inputDirBED and isPairedBED must be supplied.")
             }
              if (isPairedBAM==TRUE) {
                  tssObj@dataTypeBAM <- c("pairedEnd")
              }
              if (isPairedBAM==FALSE) {
                  tssObj@dataTypeBAM <- c("singleEnd")
              }
              
              if (isPairedBED==TRUE) {
                  tssObj@dataTypeBED <- c("pairedEnd")
              }
              if (isPairedBED==FALSE) {
                  tssObj@dataTypeBED <- c("singleEnd")
              }
              if (is.character(inputDirBAM)) {
                  tss_filesBAM <- list.files(inputDirBAM, pattern="\\.bam$",
                                          all.files=FALSE, full.names=TRUE)
                  if (length(tss_filesBAM) < 1) {
                      stop("There are no .bam files in the directory you",
                           "specified, or the directory itself does not exist.",
                           "\n Please check your input for the argument 'inputDirBAM'.")
                  }
              tssObj@fileNamesBAM <- tss_filesBAM
                  if (is.character(tssObj@dataTypeBAM)==FALSE) {
                      stop("The argument 'isPairedBAM' is empty. Please fix.")
                  }
                  if(tssObj@dataTypeBAM=="pairedEnd") {
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
                  if (tssObj@dataTypeBAM=="singleEnd") {
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
                  bam.paths <- tssObj@fileNamesBAM
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
                  #now for BED files (if present)
                  if (is.character(inputDirBED)) {
                  tss_filesBED <- list.files(inputDirBED, pattern="\\.bed$",
                                      all.files=FALSE, full.names=TRUE)
                      if (length(tss_filesBED) < 1) {
                          stop("There are no .bed files in the directory you",
                               "specified, or the directory itself does not exist.",
                               "\n Please check your input for the argument 'inputDirBED'.")
                      }
                      if (is.character(tssObj@dataTypeBED==FALSE)) {
                          stop("The argument 'isPairedBAM' is empty. Please fix.")
                      }
                  if (isPairedBED == TRUE) {
                  message("\nImporting paired-end reads ...\n")
                  tssObj@fileNamesBED <- tss_filesBED
                  bed.paths <- tssObj@fileNamesBED
                  n.beds <- length(tss_files)
                  message("\nBeginning import of ", n.beds, " bed files ...\n")
                  beds.GR <- bplapply(tss_files, import,
                                      format="bedpe",
                                      BPPARAM = MulticoreParam()
                                      )
                  tssObj@bedData <- beds.GR
                  message("Done. Alignment data from ", n.beds,
                          " bed files have been attached to the tssObject.\n")
                  message("-----------------------------------------------------\n")
              }
              if (isPairedBED==FALSE) {
              message("\nImporting single-end reads ...\n")
              tssObj@fileNamesBED <- tss_files
              bed.paths <- tssObj@fileNamesBED
              n.beds <- length(tss_files)
              message("\nBeginning import of ", n.beds, " bed files ...\n")
              beds.GR <- bplapply(bed.paths, import.bed, 
                                  BPPARAM = MulticoreParam()
                                  )
              tssObj@bedData <- beds.GR
              message("Done. Alignment data from ", n.beds,
                  " bed files have been attached to the tssObject.\n")
              message("-----------------------------------------------------\n")
          }
              }
              if (length(sampleNames)!=(length(tssObj@fileNamesBAM) + length(tssObj@fileNamesBED))) { #accounting for BAM + BED
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

