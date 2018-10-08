#' @title \strong{cDNAtoBED}
#' @description \code{cDNAtoBED} converts aligned cDNA data (in .gsq format)
#' to BED format, extracting the 5'-most base.
#'
#' @param gsqFile a path to the gsq (GeneSeqer) output file (class character)
#' @param fileName the name (class character) of the BED file to be written
#' (default is "gsqOut.bed")
#'
#' @return a BED file containing a list of the 5'-most base from each of the
#' alignments contained in the GeneSeqer (.gsq) output file is written to the
#' user's working directory.
#' 
#' @import BiocGenerics
#'
#' @examples
#' extdata.dir <- system.file("extdata", package="TSRchitect")
#' tssObjectExample <- cDNAtoBED(gsqFile=paste(extdata.dir,"AtEST.gsq",sep="/"),
#'                               fileName="testOut.bed")
#'
#' @export
#' @rdname cDNAtoBED-methods

setGeneric("cDNAtoBED",
    function(gsqFile, fileName="gsqOut.bed")
    standardGeneric("cDNAtoBED")
)

#' @rdname cDNAtoBED-methods

setMethod("cDNAtoBED",
          signature(gsqFile="character", fileName="character"),
          function(gsqFile, fileName="gsqOut.bed") {

              message("... cDNAtoBED ...")

              gsq_extractor <- system.file("extdata", "gsq_gff_extractor.pl", package="TSRchitect")
              
              arg1 <- paste("-gsq", gsqFile)
              arg2 <- paste("-out", fileName)
              cmd <- paste("perl", gsq_extractor, arg1, arg2)
              system(cmd)

              message("The file ", fileName,
                      " has been written to your working directory."
                      )

          }
          )
