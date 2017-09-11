#' getPathBed
#' @description An internal function  which is invoked using the user-level
#' function loadTSSObject to import BED data.
#' (Internal function)
#'
#' @param experimentName - a S4 object of class tssObject containing
#' information in slot tssTagData
#' @param setNum - a numeric value representing the fileName slot
#' to access
#' 
#  @keywords internal
#'
#' @return returns a single GRanges object from a BED file (internal)
#'
#' @importFrom rtracklayer import import.bed 
#' 
#' @export
#' @rdname getPathBed-methods


setGeneric("getPathBed",
    function(experimentName, setNum, is.bedpe)
    standardGeneric("getPathBed")
)

#' @rdname getPathBed-methods

setMethod("getPathBed",
          signature(experimentName="tssObject", "numeric", "logical"),

          function(experimentName, setNum, is.bedpe=FALSE) {
              message("... getPathBed ...")
              path.names <- experimentName@fileNames

              if (setNum > length(path.names)) {
                  stop("setNum must not exceed the total number of datasets.")
              }
              this.path <- path.names[setNum]

              if (is.bedpe=="TRUE") {
                  this.GR <- import(con=this.path, format="bedpe")
              }

              else {
                  this.GR <- import.bed(con=this.path)
              }
              message("-----------------------------------------------------\n")
              message(" Done.\n")
              return(this.GR)
      }
          )
