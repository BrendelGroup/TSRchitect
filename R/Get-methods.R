#' @title \strong{getTitle}
#' @description an accessor function that retrieves the contents of slot "title"
#' from a given \emph{tssObject}
#'
#' @param experimentName n S4 object of class \emph{tssObject}
#'
#' @return the contents of slot "title" are returned, which are of class
#' "character"
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
#' @rdname getTitle-methods


setGeneric("getTitle",
    function(experimentName)
    standardGeneric("getTitle")
)

#' @rdname getTitle-methods

setMethod("getTitle",
           signature(experimentName = "tssObject"),
           function (experimentName){
               my.title <- experimentName@title
               return(my.title)
	   }
)

#' @title \strong{getFileNames}
#' @description an accessor function that retrieves the contents of slot
#' "fileNames" from a given \emph{tssObject}
#'
#' @param experimentName an S4 object of class \emph{tssObject}
#'
#' @return the contents of slot "fileNames" are returned, which are
#' a vector of class "character"
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
#' @rdname getFileNames-methods

setGeneric(
name="getFileNames",
def=function(experimentName){
	standardGeneric("getFileNames")
}
)

#' @rdname getFileNames-methods

setMethod("getFileNames",
signature(experimentName = "tssObject"),
function (experimentName){
        my.files <- vector(mode="character")
    if (!is.na(experimentName@fileNamesBAM[[1]])) {
        my.files1 <- experimentName@fileNamesBAM
        my.files <- c(my.files, my.files1)
    }
    if (!is.na(experimentName@fileNamesBED[[1]])) {
        my.files1 <- experimentName@fileNamesBED
        my.files <- c(my.files, my.files1)
    }
    return(my.files)
}
)

#' @title \strong{getBamData}
#' @description an accessor function that retrieves the contents of
#' a specified slot "bamData" from a given \emph{tssObject}
#'
#' @param experimentName an S4 object of class \emph{tssObject}
#' @param slot 'numeric' a number corresponding to the slot in
#' "bamData" to be retrieved.
#'
#' @return the contents of the specified slot "bamData" are returned
#'
#' @keywords methods
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' example.bamData <- getBamData(experimentName=tssObjectExample, slot = 1)
#' example.bamData
#'
#' @export
#' @rdname getBamData-methods

setGeneric(
name="getBamData",
def=function(experimentName, slot){
	standardGeneric("getBamData")
}
)

#' @rdname getBamData-methods

setMethod("getBamData",
signature(experimentName = "tssObject", slot = "numeric"),
function (experimentName, slot){
    n.bams <- length(experimentName@bamData)
    if (slot > n.bams) {
        stop("The @bamData slot you selected does not exist.\n")
    }
    bam.data <- experimentName@bamData[[slot]]
    return(bam.data)
}
)

#' @title \strong{getTSStagData}
#' @description an accessor function that retrieves the contents of
#' a specified slot "tssTagData" from a given \emph{tssObject}
#'
#' @param experimentName an S4 object of class \emph{tssObject}
#' @param slot 'numeric' a number corresponding to the slot in
#' "tssTagData" to be retrieved.
#'
#' @return the contents of the specified slot "tssTagData" are returned
#'
#' @keywords methods
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' example.tssTagData <- getTSStagData(experimentName=tssObjectExample, slot=1)
#' example.tssTagData
#'
#' @export
#' @rdname getTSStagData-methods

setGeneric(
name="getTSStagData",
def=function(experimentName, slot){
	standardGeneric("getTSStagData")
}
)

#' @rdname getTSStagData-methods

setMethod("getTSStagData",
          signature(experimentName = "tssObject", slot = "numeric"),
          function (experimentName, slot){

              n.tss <- length(experimentName@tssTagData)
              if (slot > n.tss) {
                  stop("The @tssTagData slot you selected does not exist.\n")
              }
              tss.data <- experimentName@tssTagData[[slot]]
              return(tss.data)
}
)

#' @title \strong{getTSScountData}
#' @description an accessor function that retrieves the contents of
#' a specified slot "tssCountData"/"tssCountDataMerged" from a
#' given \emph{tssObject}
#'
#' @param experimentName an S4 object of class \emph{tssObject}
#' @param slotType 'character' which data type is to be selected.
#' Either "replicates" (tssCountData) or "merged" (tssCountDataMerged)
#' @param slot 'numeric' a number corresponding to the slot in
#' "tssCountData"/"tssCountDataMerged" to be retrieved.
#'
#' @return the contents of the selected slot (either "tssCountData" or
#' "tssCountDataMerged" are returned)
#'
#' @keywords methods
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' ex.tssCountData <- getTSScountData(experimentName=tssObjectExample,
#' slotType="replicates", slot = 1)
#' ex.tssCountData
#'
#' @export
#' @rdname getTSScountData-methods

setGeneric(
name="getTSScountData",
def=function(experimentName, slotType, slot){
	standardGeneric("getTSScountData")
}
)

#' @rdname getTSScountData-methods

setMethod("getTSScountData",
          signature(experimentName = "tssObject", slotType = "character",
                    slot = "numeric"),
          function (experimentName, slotType = c("replicates", "merged"),
                    slot) {
              fileType <- match.arg(slotType, several.ok=FALSE)
              if (slotType=="replicates") {
                  n.tss.counts  <- length(experimentName@tssCountData)
                  if (slot > n.tss.counts) {
                      stop("The @tssTagData slot you selected",
                           " does not exist.\n")
                  }
                  my.counts.data <- experimentName@tssCountData[[slot]]
              }

              if (slotType=="merged") {
                  n.tss.counts  <- length(experimentName@tssCountDataMerged)
                  if (slot > n.tss.counts) {
                      stop("The @tssTagDataMerged slot you selected",
                           "does not exist.\n")
                  }
                  my.counts.data <- experimentName@tssCountDataMerged[[slot]]
              }
              return(my.counts.data)
          }
)

#' @title \strong{getTSRdata}
#' @description an accessor function that retrieves the contents of
#' a specified slot "tsrData"/"tsrDataMerged" from a
#' given \emph{tssObject}
#'
#' @param experimentName an S4 object of class \emph{tssObject}
#' @param slotType 'character' which data type is to be selected.
#' Either "replicates" (tsrCountData) or "merged" (tsrCountDataMerged)
#' @param slot 'numeric' a number corresponding to the slot in
#' "tsrData"/"tsrDataMerged" to be retrieved.
#'
#' @return the contents of the selected slot (either "tsrData" or
#' "tsrDataMerged" are returned)
#'
#' @keywords methods
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' ex.tsrData <- getTSRdata(experimentName=tssObjectExample,
#' slotType="replicates", slot = 1)
#' ex.tsrData
#'
#' @export
#' @rdname getTSRdata-methods

setGeneric(
name="getTSRdata",
def=function(experimentName, slotType, slot){
	standardGeneric("getTSRdata")
    }
)

#' @rdname getTSRdata-methods

setMethod("getTSRdata",
          signature(experimentName = "tssObject", slotType = "character",
                    slot = "numeric"),
          function (experimentName, slotType = c("replicates", "merged"),
                    slot) {
                    fileType <- match.arg(slotType, several.ok=FALSE)
                    if (slotType=="replicates") {
                        n.tsrs  <- length(experimentName@tsrData)
                        if (slot > n.tsrs) {
                            stop("The @tsrData slot you selected",
                                 "does not exist.\n")
                     }
                  my.tsrs <- experimentName@tsrData[[slot]]
                  }

              if (slotType=="merged") {
                  n.tsrs  <- length(experimentName@tsrDataMerged)
                  if (slot > n.tsrs) {
                      stop("The @tsrDataMerged slot you selected",
                           "does not exist.\n")
                  }
                  my.tsrs <- experimentName@tsrDataMerged[[slot]]
              }
              return(my.tsrs)
   }
)
