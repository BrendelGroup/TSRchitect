#' Computes TSS positions from all from a tssExp object
#' @param expName an S4 object of class tssExp that contains information about the experiment.
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom Rsamtools scanBamFlag ScanBamParam BamViews
#' @return bam data from the tss profiling experiments assigned to your tssExp object (expName)
#' @export 

setGeneric(
    name="importBam",
    def=function(expName) {
        standardGeneric("importBam")
    }
    )

setMethod("importBam",
          signature(expName="tssExp"),
          function(expName) {
              expName.chr <- deparse(substitute(expName))
              exp.type <- expName@dataType

             if(exp.type=="pairedEnd") {
                  scanBamFlag(isPaired=TRUE, isProperPair=TRUE, isUnmappedQuery=FALSE) -> bamFlags
                  cat("\nPaired-end TSS data specified.")
                  c("rname","strand","pos","qwidth", "mapq", "isize") -> myFields
              }
                                     
              else {
                  scanBamFlag(isPaired=FALSE, isProperPair=FALSE, isUnmappedQuery=FALSE) -> bamFlags
                  cat("\nSingle-end TSS data specified.\n")
                  c("rname","strand","pos", "qwidth", "mapq") -> myFields
              }

              my.param <- ScanBamParam(flag=bamFlags, what=myFields)
              bam.paths <- expName@fileNames
              bv_obj <- BamViews(bam.paths)
              bv_files <- dimnames(bv_obj)[[2]]
              n.bams <- length(bv_files)
              message("\nBegan import of", n.bams, "bam files.\n")
              bams.GA <- bplapply(bam.paths, readGAlignments, BPPARAM = MulticoreParam(),param=my.param)
              expName@bamData <- bams.GA
              message("\nImport complete!\n")
              message("\nAlignment data from", n.bams, "bams have been attached to your tssExp object.\n")
              assign(expName.chr, expName, parent.frame()) 
          }
          )
