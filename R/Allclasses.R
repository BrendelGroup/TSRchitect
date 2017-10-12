setClass(Class="tssObject",
         representation(
             title = "character",
             inputDir = "character",
             fileNamesBAM = "character",
             fileNamesBED = "character",
             dataTypeBAM = "character",
             dataTypeBED = "character",
             sampleNames = "character",
             replicateIDs = "numeric",
             annotation = "GRanges",
             bamData = "list",
             bedData = "list",
             tssTagData = "GRangesList",
             tssCountData = "list",
             tssCountDataMerged = "list",
             countsData = "data.frame",
             tsrData = "list",
             tsrDataMerged = "list",
             tsrCounts = "data.frame"
             ),
         prototype(
             title = NA_character_,
             inputDir = NA_character_,
             fileNamesBAM = NA_character_,
             fileNamesBED = NA_character_,
             dataTypeBAM = NA_character_,
             dataTypeBED = NA_character_,
             sampleNames = NA_character_,
             replicateIDs = vector("integer"),
             annotation = GRanges(),
             bamData = list(),
             bedData = list(),
             tssTagData = GRangesList(),
             tssCountData = list(),
             tssCountDataMerged = list(),
             countsData = data.frame(),
             tsrData = list(),
             tsrDataMerged = list(),
             tsrCounts = data.frame()
             ),
         )

#' @title \strong{tssObject}
#' @description S4 constructor function for \emph{tssObject}
#'
#' @param title 'character' A short descriptive title for the
#' experiment. Is set to NA by default.
#' @param bamData 'list' the name of a list of \linkS4class{GAlignments}
#' objects (originating from .bam files) in the workspace. 
#' Set to NA by default.
#' @param bedData 'list' the name of a list of \linkS4class{GRanges} or
#' \linkS4class{Pairs} objects (originating from .bed files) in the
#' workspace. Set to NA by default.
#'
#' @return a new \emph{tssObject} is returned to the user's workspace.
#'
#' @importClassesFrom S4Vectors Pairs
#'
#' @examples
#' new.tssObj <- tssObject(title="Example")
#'
#' @export

tssObject <- function(title=NA, bamData=NA, bedData=NA) {

    new.tssObj <- new("tssObject")

    if (!(is.na(title))) {
        new.tssObj@title <- title
    }

    if (!(is.na(bamData))) {
        new.tssObj@bamData <- bamData
    }

    if (!(is.na(bedData))) {
        new.tssObj@bedData <- bedData
    }

    return(new.tssObj)
}
