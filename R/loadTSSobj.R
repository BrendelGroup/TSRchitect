#' @title \strong{loadTSSobj}
#' @description \code{loadTSSobj} processes alignment files in .bam or .bed
#' formats from the local directory supplied.
#'
#' @param experimentTitle a descriptive title for the experiment (character).
#' @param inputDir path to the directory containing the alignment files (in
#' either .bam or .bed formats) (character). 
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bam or .bed will be imported with this function.
#' @param isPairedBAM if the input is in BAM format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE) (logical)
#' @param isPairedBED if the input is in BED format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE). Set to FALSE by default. (logical)
#' Note: if TRUE, the input data must be in bedpe format, as described here:
#' http://bedtools.readthedocs.io/en/latest/content/general-usage.html
#' @param sampleNames unique labels of class character for each TSS sample
#' within the experiment (character).
#' @param replicateIDs identifiers indicating which samples are biological
#' replicates. Note that \code{loadTSSobj} imports
#' alignment data in ascending alphanumeric order, so the arguments to
#' replicateIDs must be arranged in this order also so that they directly
#' correspond to the intended file (numeric).
#'
#' @return \emph{loadTSSobj} fills the slot \emph{bamData} and/or \emph{bedData}
#' on the returned \emph{tssObject} with \linkS4class{GAlignments} objects
#' (for .bam files), or \linkS4class{GRanges} objects (for .bed files).
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#' @importFrom methods new
#' @importFrom rtracklayer import import.bed
#'
#' @examples
#' extdata.dir <- system.file("extdata", package="TSRchitect")
#' test.Obj <- loadTSSobj(experimentTitle="Code example", inputDir=extdata.dir,
#' isPairedBAM=TRUE, sampleNames=c("sample1-rep1", "sample1-rep2",
#' "sample2-rep1","sample2-rep2"), replicateIDs=c(1,1,2,2))
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' @note All files found in \emph{inputDir} will be
#' retrieved and written in ascending alphanumeric order to the
#' \emph{@fileNamesBAM} and/or \emph{@fileNamesBED} slot(s) on
#' the \emph{tssObject} that is created.
#'
#' @export
#' @rdname loadTSSobj-methods


setGeneric("loadTSSobj",
           function(experimentTitle, inputDir, 
                    isPairedBAM=FALSE, isPairedBED=FALSE, 
                    sampleNames,
                    replicateIDs)
    standardGeneric("loadTSSobj")
)

#' @rdname loadTSSobj-methods

setMethod("loadTSSobj", #both BAM and BED files
          signature(experimentTitle="character", inputDir="character",
                    isPairedBAM="ANY", isPairedBED="ANY",
                    sampleNames="character",
                    replicateIDs="numeric"),
          function(experimentTitle, inputDir, isPairedBAM=FALSE,
                   isPairedBED=FALSE,
                   sampleNames,
                   replicateIDs) {
              
              message("... loadTSSobj ...")
              tssObj <- new("tssObject")

              tssObj@title <- experimentTitle

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
             tss_filesBAM <- list.files(inputDir, pattern="\\.bam$",
                                        all.files=FALSE, full.names=TRUE)
             tssObj@fileNamesBAM <- tss_filesBAM
             if (length(tss_filesBAM) > 0) {
                  if (is.character(tssObj@dataTypeBAM)==FALSE) {
                      stop("The argument 'isPairedBAM' is empty. Please fix.")
                  }
                  if (tssObj@dataTypeBAM=="pairedEnd") {
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
              tss_filesBED <- list.files(inputDir, pattern="\\.bed$",
                                         all.files=FALSE, full.names=TRUE)
              if (length(tss_filesBED) > 0) {
                  if (isPairedBED == TRUE) {
                  message("\nImporting paired-end reads ...\n")
                  tssObj@fileNamesBED <- tss_filesBED
                  bed.paths <- tssObj@fileNamesBED
                  n.beds <- length(tss_filesBED)
                  message("\nBeginning import of ", n.beds, " bed files ...\n")
                  beds.GR <- bplapply(tss_filesBED, import,
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
              tssObj@fileNamesBED <- tss_filesBED
              bed.paths <- tssObj@fileNamesBED
              n.beds <- length(tss_filesBED)
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
#              if (length(sampleNames)!=(length(tssObj@fileNamesBAM) + length(tssObj@fileNamesBED))) { #accounting for BAM + BED
#                  print(length(tssObj@fileNamesBED))
#                  stop("\nNumber of sampleNames must be equal to",
#                       " number of input files.")
#              }
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

