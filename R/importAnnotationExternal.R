#' @title \strong{importAnnotationExternal}
#' @description \code{importAnnotationExternal} imports an
#' annotation from an external file and attaches it to the
#' \emph{tssObject}
#'
#' @param experimentName - an S4 object of class \emph{tssObject}
#' that contains information about the experiment
#' @param fileType - the format of the annotation file to be imported.
#' Must be one of: "bed", "gff" or "gff3".
#' @param annotFile - a path (full or relative) to the annotation
#' file to be imported.
#'
#' @return fills the slot \emph{@@annotation} in the returned
#' \emph{tssObject} with a \linkS4class{GRanges} object contining a
#' parsed annotation file of the selected type.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.bed import.gff import.gff3
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))
#' extdata.dir <- system.file("extdata", package="TSRchitect")
#' annotation <- dir(extdata.dir, pattern="\\.gff3$", full.names=TRUE)
#' tssObjectExample <- importAnnotationExternal(experimentName=tssObjectExample,
#' fileType="gff3", annotFile=annotation)
#'
#' @note \code{importAnnotationExternal} makes use of three functions from the
#' \emph{rtracklayer} package: \code{\link[rtracklayer]{import.bed}},
#' \code{\link[rtracklayer]{import.gff}}, and
#' \code{\link[rtracklayer]{import.gff3}}
#' @note An example similar to the one provided can be found in
#' \emph{Example 2} from the vignette (/inst/doc/TSRchitect.Rmd)
#' @note To import directly from the AnnotationHub client, please
#' refer to \code{importAnnotationHub}
#'
#' @export
#' @rdname importAnnotationExternal-methods


setGeneric("importAnnotationExternal",
    function(experimentName, fileType, annotFile)
    standardGeneric("importAnnotationExternal")
)

#' @rdname importAnnotationExternal-methods

setMethod("importAnnotationExternal",
          signature(experimentName="tssObject", fileType="character",
                    annotFile="character"),
          function(experimentName, fileType=c("bed","gff","gff3"),
                   annotFile) {

              message("... importAnnotationExternal ...")
              fileType <- match.arg(fileType, several.ok=FALSE)
              if (fileType=="bed") {
                  experimentName@annotation <- import.bed(annotFile)
                  }
              if (fileType=="gff") {
                  experimentName@annotation <- import.gff(annotFile)
                  }
              if (fileType=="gff3") {
                  experimentName@annotation <- import.gff3(annotFile)
                  }
              message("Done. Annotation data have been attached to",
                  " the tssObject.\n")
              message("-----------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
