setClass(Class="tssObject",
         representation(
             title = "character",
             fileNames = "character",
             dataType = "character",
             sampleNames = "character",
             replicateIDs = "numeric",
             annotation = "GRanges",
             bamData = "list",
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
             fileNames = NA_character_,
             dataType = NA_character_,
             sampleNames = NA_character_,
             replicateIDs = vector("integer"),
             annotation = GRanges(),
             bamData = list(),
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
#' objects in the workspace. Set to NA by default.
#'
#' @return a new \emph{tssObject} is returned to the user's workspace.
#' 
#' @examples
#' new.tssObj <- tssObject(title="Example")
#'
#' @export

tssObject <- function(title=NA, bamData=NA) {

    new.tssObj <- new("tssObject")

    if (!(is.na(title))) {
        new.tssObj@title <- title
    }

    if (!(is.na(bamData))) {
        new.tssObj@bamData <- bamData
    }

    message("\nA new tssObject has been created in your workspace.\n")
    return(new.tssObj)
}
                      
