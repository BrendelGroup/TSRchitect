#' exportClasses(TSRchitect)

setGeneric(
name="getTitle",
def=function(expName){
	standardGeneric("getTitle")
}
)

setMethod("getTitle",
signature(expName = "tssExp"),
function (expName){
	expName@title
}
)

setGeneric(
name="getFileNames",
def=function(expName){
	standardGeneric("getFileNames")
}
)

setMethod("getFileNames",
signature(expName = "tssExp"),
function (expName){
    expName@fileNames
}
)

setGeneric(
name="getBamData",
def=function(expName){
	standardGeneric("getBamData")
}
)

setMethod("getBamData",
signature(expName = "tssExp"),
function (expName){
	expName@bamData
}
)

setGeneric(
name="getTSSs",
def=function(expName){
	standardGeneric("getTSSs")
}
)

setMethod("getTSSs",
signature(expName = "tssExp"),
function (expName){
	expName@tssData
}
)

setGeneric(
name="getExpr",
def=function(expName){
	standardGeneric("getExpr")
}
)

setMethod("getExpr",
signature(expName = "tssExp"),
function (expName){
	expName@expData
}
)

setGeneric(
name="getTSRs",
def=function(expName){
	standardGeneric("getTSRs")
    }
)

setMethod("getTSRs",
signature(expName = "tssExp"),
function (expName){
	expName@tsrData
}
)
