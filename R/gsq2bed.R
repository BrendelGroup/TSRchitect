#' @title \strong{gsq2bed}
#' @description \code{gsq2bed} converts aligned cDNA data (in .gsq format)
#' to BED format, extracting the 5'-most base.
#'
#' @param gsqFile a path to the gsq (GeneSeqer) output file (class character)
#' @param outfile the name (class character) of the BED file to be written
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
#' tssObjectExample <- gsq2bed(gsqFile=paste(extdata.dir,"AtEST.gsq",sep="/"),
#'                             outfile="")
#'
#' @export
#' @rdname gsq2bed-methods

setGeneric("gsq2bed",
    function(gsqFile, outfile)
    standardGeneric("gsq2bed")
)

#' @rdname gsq2bed-methods

setMethod("gsq2bed",
          signature(gsqFile="character", outfile="character"),
          function(gsqFile, outfile) {

              message("... gsq2bed ...")

	      if (missing(outfile)) {outfile=paste(gsqFile,"bed",sep=".") }
              perl_script <- system.file("extdata/Perl", "gsq2bed.pl",
                                          package="TSRchitect")
              arg1 <- paste("-gsq", gsqFile)
              arg2 <- paste("-out", outfile)
              cmd  <- paste("perl", perl_script, arg1, arg2)
              system(cmd)

              message("The file ", outfile,
                      " has been generated in your working directory."
                      )

          }
          )
