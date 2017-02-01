#' @title \strong{getTitle}
#' @description an accessor function that retrieves the contents of slot "title"
#' from a given \emph{tssObject}
#'
#' @param experimentName n S4 object of class \emph{tssObject}
#' 
#' @return the contents of slot "title" are returned
#'
#' @keywords methods
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' example.title <- getTitle(experimentName=tssObjectExample)
#' example.title
#'
#' @export

setGeneric(
name="getTitle",
def=function(experimentName){
	standardGeneric("getTitle")
}
)

setMethod("getTitle",
signature(experimentName = "tssObject"),
function (experimentName){
	experimentName@title
}
)

#' @title \strong{getFileNames}
#' @description an accessor function that retrieves the contents of slot "fileNames"
#' from a given \emph{tssObject}
#'
#' @param experimentName n S4 object of class \emph{tssObject}
#' 
#' @return the contents of slot "fileNames" are returned
#'
#' @keywords methods
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' example.file.names <- getFileNames(experimentName=tssObjectExample)
#' example.file.names
#'
#' @export

setGeneric(
name="getFileNames",
def=function(experimentName){
	standardGeneric("getFileNames")
}
)

setMethod("getFileNames",
signature(experimentName = "tssObject"),
function (experimentName){
    experimentName@fileNames
}
)

setGeneric(
name="getBamData",
def=function(experimentName){
	standardGeneric("getBamData")
}
)

setMethod("getBamData",
signature(experimentName = "tssObject"),
function (experimentName){
	experimentName@bamData
}
)

setGeneric(
name="getTSSs",
def=function(experimentName){
	standardGeneric("getTSSs")
}
)

setMethod("getTSSs",
signature(experimentName = "tssObject"),
function (experimentName){
	experimentName@tssTagData
}
)

setGeneric(
name="getExpr",
def=function(experimentName){
	standardGeneric("getExpr")
}
)

setMethod("getExpr",
signature(experimentName = "tssObject"),
function (experimentName){
	experimentName@tssCountData
}
)

setGeneric(
name="getTSRs",
def=function(experimentName){
	standardGeneric("getTSRs")
    }
)

setMethod("getTSRs",
signature(experimentName = "tssObject"),
function (experimentName){
	experimentName@tsrData
}
)
