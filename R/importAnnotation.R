#' @title \emph{importAnnotation}
#'
#' @description \emph{importAnnotation} processes BAM files as specified by \emph{initializeExp}.
#'
#' @param experimentName - an S4 object of class tssObject that contains information about the experiment
#'
#' @param fileType - the format of the annotation file to be imported. Must be one of: "bed", "gff" or "gff3".
#'
#' @param annotFile - a path (full or relative) to the annotation file to be imported. 
#' 
#' @return \emph{impartAnnotation} fills the slot experimentName@geneAnnot in the tssObject \emph{experimentName} with
#'         a GRanges object contining the annotion data
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.bed import.gff import.gff3
#' @export

setGeneric(
    name="importAnnotation",
    def=function(experimentName, fileType, annotFile) {
        standardGeneric("importAnnotation")
    }
    )

setMethod("importAnnotation",
          signature(experimentName="tssObject", fileType="character", annotFile="character"),
          function(experimentName, fileType=c("bed","gff3"), annotFile) {
              experimentName.seq <- deparse(substitute(experimentName))
              message("... importAnnotation ...")
              fileType <- match.arg(fileType, c("bed","gff", "gff3"), several.ok=FALSE)
              if (fileType=="bed") {
                  import.bed(annotFile) -> my.annot
                  }
              if (fileType=="gff") {
                  import.gff(annotFile) -> my.annot
                  }
              if (fileType=="gff3") {
                  import.gff3(annotFile) -> my.annot
                  }
              experimentName@geneAnnot <- my.annot
              cat("Done. Annotation data has been attached to tssObject\nobject \"", experimentName.seq, "\".\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(experimentName.seq, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
