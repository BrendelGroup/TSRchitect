#' @title \strong{importAnnotationHub}
#' @description imports an annotation directly from AnnotationHub
#' and attaches it to the \emph{tssObject}
#'
#' @param experimentName - an S4 object of class \emph{tssObject}
#' that contains information about the experiment
#' @param provider - 'character' the source of the annotation in AnnotationHub
#' (e.g. 'gencode')
#' @param annotType - 'character' the format of the annotationHub annotation
#' to be imported. Using 'gff' will find annotations in both gff or gff3 formats
#' @param species - 'character' the species identifier in AnnotationHub
#' (e.g. 'human')
#' @param annotID - 'character' the AnnotationHub identifier to be retrieved
#'
#' @return fills the slot \emph{@@annotation} in the returned \emph{tssObject}
#' with an AnnotationHub record. The record retrieved must be an object
#' of class \linkS4class{GRanges}.
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd)
#' @note Please consult the available records in AnnotationHub
#' beforehand using \code{\link[AnnotationHub]{AnnotationHub-class}}
#'
#' @export
#' @rdname importAnnotationHub-methods


setGeneric("importAnnotationHub",
    function(experimentName, provider, annotType, species, annotID)
    standardGeneric("importAnnotationHub")
)

#' @rdname importAnnotationHub-methods

setMethod("importAnnotationHub",
          signature(experimentName="tssObject", provider="character",
                    annotType="character", species="character",
                    annotID="character"),
          function(experimentName, provider, annotType, species,
                   annotID) {

              message("... importAnnotationHub ...")
              hub <- AnnotationHub()
              query(hub, c(provider, annotType, species))
              message("\nRetrieving selected record from AnnotationHub record ",
                      annotID, ".\n")
              annot.object <- hub[[annotID]]
              if (class(annot.object) != "GRanges") {
                  stop("\nThe selected annotation record is not",
                       " a GRanges object. Please select another.")
              }
              else {
                  experimentName@annotation <- annot.object
              }
              message("Done. Annotation data have been attached to",
                  " the tssObject.\n")
              message("---------------------------------------------------------\n")
              message(" Done.\n")
              return(experimentName)
          }
          )
