#' @title \strong{importAnnotation}
#' @description \code{importAnnotation} imports an annotation from an external file and attaches it to \emph{experimentName}.
#' @param experimentName - an S4 object of class \emph{tssObject} that contains information about the experiment
#' @param fileType - the format of the annotation file to be imported. Must be one of: "bed", "gff" or "gff3".
#' @param annotFile - a path (full or relative) to the annotation file to be imported. 
#' @return \emph{importAnnotation} fills the slot \emph{@@annotation} in the \emph{tssObject} with
#'         a \linkS4class{GRanges} object contining parsed annotation file of the selected type.
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.bed import.gff import.gff3
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' load(system.file("extdata", "genocode.v19.annotation.gff3", package="TSRchitect"))
#' importAnnotation(experimentName=tssObjectExample, fileType ="gff3", annotFile ="data/gencode.v19.annotation.gff3")
#' @note An example similar to the one provided can be found in \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd)
#' @note \code{importAnnotation} makes use of three functions from the \emph{rtracklayer} package: \code{\link[rtracklayer]{import.bed}}, \code{\link[rtracklayer]{import.gff}}, and \code{\link[rtracklayer]{import.gff3}}
#' @export

setGeneric(
    name="importAnnotation",
    def=function(experimentName, fileType, annotFile) {
        standardGeneric("importAnnotation")
    }
    )

setMethod("importAnnotation",
          signature(experimentName="tssObject", fileType="character", annotFile="character"),
          function(experimentName, fileType=c("bed","gff","gff3"), annotFile) {
              object.name <- deparse(substitute(experimentName))
              message("... importAnnotation ...")
              fileType <- match.arg(fileType, c("bed","gff", "gff3"), several.ok=FALSE)
              if (fileType=="bed") {
                  import.bed(annotFile) -> experimentName@annotation
                  }
              if (fileType=="gff") {
                  import.gff(annotFile) -> experimentName@annotation
                  }
              if (fileType=="gff3") {
                  import.gff3(annotFile) -> experimentName@annotation
                  }
              cat("Done. Annotation data has been attached to tssObject\nobject \"", object.name, "\".\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
