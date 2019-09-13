#' @title \strong{loadTSSobj}
#' @description \code{loadTSSobj} processes alignment files in .bam or .bed
#' formats from the local directory supplied.
#'
#' @param experimentTitle a descriptive title for the experiment (character).
#' @param inputDir path to the directory containing the alignment files (in
#' either .bam or .bed formats) (character). 
#' Note that all the paths to all files in \emph{inputDir} with the extension
#' .bam or .bed will be imported with this function.
#' @param n.cores the number of cores to be used for this job. (numeric)
#' @param isPairedBAM if the input is in BAM format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE). Set to TRUE by default. (logical)
#' @param isPairedBED if the input is in BED format, specifies whether the
#'  TSS profiling experiment is paired-end (if TRUE) or single-end
#'  (if FALSE). Set to TRUE by default. (logical)
#' Note: if TRUE, the input data must be in bedpe format, as described here:
#' http://bedtools.readthedocs.io/en/latest/content/general-usage.html
#' @param sampleSheet file providing TSS sample information; if provided,
#'  sampleNames and replicateIDs are ignored; input format is tab-delimited
#'  3-column rows with sample name, replicate ID, and sample data file name;
#'  the file has to include column headers "SAMPLE	ReplicateID	FILE"
#'  and can be either tab-delimited or an Excel spreadsheet (file extension
#'  ".xls" or "xlsx" (character)
#' @param sampleNames unique labels of class character for each TSS sample
#' within the experiment (character).
#' @param replicateIDs identifiers indicating which samples are biological
#' replicates. Note that \code{loadTSSobj} imports alignment data in ascending
#' alphanumeric order, so the arguments to replicateIDs must be arranged in this
#' order also so that they directly correspond to the intended file (numeric).
#'
#' @return \emph{loadTSSobj} fills the slot \emph{bamDataFirstRead} and/or \emph{bedData}
#' on the returned \emph{tssObject} with \linkS4class{GAlignments} objects
#' (for .bam files), or \linkS4class{GRanges} objects (for .bed files).
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#' @importFrom methods new
#' @importFrom readxl read_excel
#' @importFrom rtracklayer import import.bed
#' @importFrom tools file_ext
#' @importFrom utils read.table
#'
#' @examples
#' extdata.dir <- system.file("extdata/bamFiles", package="TSRchitect")
#' test.Obj <- loadTSSobj(experimentTitle="Code example", inputDir=extdata.dir,
#' n.cores=2, isPairedBAM=TRUE, sampleNames=c("sample1-rep1", "sample1-rep2",
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
           function(experimentTitle, inputDir, ...)
           standardGeneric("loadTSSobj")
)

#' @rdname loadTSSobj-methods

setMethod("loadTSSobj",
          signature(experimentTitle="character", inputDir="character"),
          function(experimentTitle, inputDir, n.cores=1,
                   isPairedBAM=TRUE, isPairedBED=TRUE,
                   sampleSheet=NA, sampleNames=NA, replicateIDs=NA) {
              
              message("... loadTSSobj ...")
              tssObj <- new("tssObject")

              tssObj@title <- experimentTitle
              tss_filesBAM <- vector(mode="character",length=0)
              tss_filesBED <- vector(mode="character",length=0)

              if(!missing(sampleSheet) & is.character(sampleSheet)) {
                sampleSheet <- file.path(inputDir,sampleSheet)
                ext <- file_ext(sampleSheet)
                if (ext %in% c("xls","xlsx")) {
                  samples <- read_excel(sampleSheet)
                } else {
                  samples <- read.table(sampleSheet,sep='	',
                                        stringsAsFactors=F,header=T, 
                                        comment.char="")
                }
                sampleNames <- samples$SAMPLE
                replicateIDs <- samples$ReplicateID
                for (file in samples$FILE) {
                  ext <- file_ext(file)
                  if (ext %in% c("bam","Bam")) {
                   tss_filesBAM <- c(tss_filesBAM,file)
                  } else if (ext %in% c("bed","Bed")) {
                   tss_filesBED <- c(tss_filesBED,file)
                  } else {
                    message("No .bam nor .bed files specified in ",sampleSheet);
                    message("Please check your input file.")
                  }
                }
              }

              if (is.na(sampleNames)  || !is.character(sampleNames) ||
                  is.na(replicateIDs) || !is.numeric(replicateIDs)    ) {
                stop("You must specify both sample names (as \"character\")",
                     " and replicate IDs (as \"numeric\").  Please check.")
              }
              if (!missing(n.cores) & n.cores > 1) {
                  BiocParallel::register(MulticoreParam(workers=n.cores),
                                                        default=TRUE)
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

              if (length(tss_filesBAM) == 0) {
                  tss_filesBAM <- list.files(inputDir, pattern="\\.bam$",
                                             all.files=FALSE, full.names=TRUE)
              }
              tssObj@fileNamesBAM <- tss_filesBAM
              if (length(tss_filesBAM) > 0) {
                if (is.character(tssObj@dataTypeBAM)==FALSE) {
                  stop("The argument 'isPairedBAM' is empty. Please fix.")
                }
                if (tssObj@dataTypeBAM=="singleEnd") {
                  message("\nImporting single-end reads ...\n")
                  bamFlags <- scanBamFlag(isPaired=FALSE,
                                          isUnmappedQuery=FALSE,
                                          isSecondaryAlignment=FALSE)
                  myFields <- c("rname","flag","pos", "qwidth","mapq","cigar")
                }

                if (tssObj@dataTypeBAM=="pairedEnd") {
                  message("\nImporting paired-end reads (first reads) ...\n")
                  bamFlags <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
                                          isFirstMateRead=TRUE, hasUnmappedMate=FALSE,
                                          isUnmappedQuery=FALSE,
                                          isSecondaryAlignment=FALSE)
                  myFields <- c("rname","flag","pos","qwidth","mapq",
                                "cigar","isize")
                }
                my.param <- ScanBamParam(flag=bamFlags, what=myFields)
                bam.paths <- tssObj@fileNamesBAM
                bv_obj <- BamViews(bam.paths)
                bv_files <- dimnames(bv_obj)[[2]]
                n.bams <- length(bv_files)
                message("\nBeginning import of ", n.bams, " bam files ...\n")
                if (!missing(n.cores) & n.cores > 1) {
                    bams.GA <- bplapply(bam.paths, readGAlignments,
                                        BPPARAM = MulticoreParam(), param=my.param)
                } else {
                    bams.GA <- lapply(bam.paths, readGAlignments, param=my.param)
                }
                tssObj@bamDataFirstRead <- bams.GA

                if (tssObj@dataTypeBAM=="pairedEnd") {
                  message("\nImporting paired-end reads (last reads) ...\n")
                  bamFlags <- scanBamFlag(isPaired=TRUE, isProperPair=TRUE,
                                          isFirstMateRead=FALSE, hasUnmappedMate=FALSE,
                                          isUnmappedQuery=FALSE,
                                          isSecondaryAlignment=FALSE)
                  my.param <- ScanBamParam(flag=bamFlags, what=myFields)
                  if (!missing(n.cores) & n.cores > 1) {
                      bams.GA <- bplapply(bam.paths, readGAlignments,
                                           BPPARAM = MulticoreParam(), param=my.param)
                  } else {
                      bams.GA <- lapply(bam.paths, readGAlignments, param=my.param)
                  }
                  tssObj@bamDataLastRead <- bams.GA
                }

                message("Done. Alignment data from ", n.bams,
                        " bam files have been attached to the tssObject.\n")
                message("-----------------------------------------------------\n")
              }

              if (length(tss_filesBED) == 0) {
                  tss_filesBED <- list.files(inputDir, pattern="\\.bed$",
                                             all.files=FALSE, full.names=TRUE)
              }
              if (length(tss_filesBED) > 0) {
                if (isPairedBED == TRUE) {
                  message("\nImporting paired-end reads ...\n")
                  tssObj@fileNamesBED <- tss_filesBED
                  n.beds <- length(tss_filesBED)
                  message("\nBeginning import of ", n.beds, " bed files ...\n")
                  if (missing(n.cores) | n.cores > 1) {
                      beds.GR <- bplapply(tss_filesBED, import, format="bedpe",
                                          BPPARAM = MulticoreParam() )
                  } else {
                      beds.GR <- lapply(tss_filesBED, import, format="bedpe")
                  }
                  tssObj@bedData <- beds.GR
                  message("Done. Alignment data from ", n.beds,
                          " bed files have been attached to the tssObject.\n")
                  message("-----------------------------------------------------\n")
                }
                if (isPairedBED==FALSE) {
                  message("\nImporting single-end reads ...\n")
                  tssObj@fileNamesBED <- tss_filesBED
                  n.beds <- length(tss_filesBED)
                  message("\nBeginning import of ", n.beds, " bed files ...\n")
                  if (!missing(n.cores) & n.cores > 1) {
                      beds.GR <- bplapply(tss_filesBED, import.bed, 
                                          BPPARAM = MulticoreParam() )
                  } else {
                      beds.GR <- lapply(tss_filesBED, import.bed)
                  }
                  tssObj@bedData <- beds.GR
                  message("Done. Alignment data from ", n.beds,
                          " bed files have been attached to the tssObject.\n")
                  message("-----------------------------------------------------\n")
                }
              }
               if (length(sampleNames)!=(length(tss_filesBAM) +
                                         length(tss_filesBED)  )) {
                  stop("\nNumber of sampleNames must be equal to",
                        " number of input files.\nPlease check.")
              }
              if (length(sampleNames)!=(length(replicateIDs))) {
                stop("\nsampleNames and replicateIDs must have equal lengths.",
                     "  Please check.")
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
