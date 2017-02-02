#' @title \strong{importAnnotationHub}
#' @description imports an annotation directly from AnnotationHub
#' and attaches it to the \emph{tssObject}
#'
#' @param experimentName - an S4 object of class \emph{tssObject}
#' that contains information about the experiment
#' @param provider - 'character' the source of the annotation in AnnotationHub
#' (e.g. 'Gencode')
#' @param annotType - 'character' the format of the annotationHub annotation
#' to be imported. Set to 'gff' by default.
#' @param species - 'character' the species identifier in AnnotationHub
#' (e.g. 'human')
#' @param annotID - 'character' the identifier corresponding to the 
#' 
#' @return fills the slot \emph{@@annotation} in the \emph{tssObject}
#' with an AnnotationHub record. The record retrieved must be an object
#' of class \linkS4class{GRanges}.
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @note An example similar to the one provided can be found in
#' \emph{Example 1} from the vignette (/inst/doc/TSRchitect.Rmd)
#' @note Please consult the available records in AnnotationHub
#' beforehand using \code{\link[AnnotationHub]{AnnoationHub}}
#'
#' @export

setGeneric(
    name="importAnnotationHub",
    def=function(experimentName, provider, annotType, species, annotID) {
        standardGeneric("importAnnotationHub")
    }
    )

setMethod("importAnnotationHub",
          signature(experimentName="tssObject", provider="character",
                    annotType="character", species="character",
                    annotID="character"),
          function(experimentName, provider, annotType, species,
                   annotID) {
              object.name <- deparse(substitute(experimentName))
              message("... importAnnotationHub ...")
              AnnotationHub() -> hub
              query(hub, c(provider, annotType, species))
              message("\nRetrieving selected record from AnnotationHub record",
                      annotID, ".\n")
              annot.object <- hub[[annotID]]
              if (class(annot.object) != "GRanges") {
                  stop("\nThe selected annotation record is not",
                       " a GRanges object. Please select another.")
              }
              else {
                  annot.object -> experimentName@annotation
              }
              cat("Done. Annotation data has been attached to",
                  "tssObject\nobject \"", object.name, "\".\n")
              cat("-------------------------------------------------------\n")
              assign(object.name, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
