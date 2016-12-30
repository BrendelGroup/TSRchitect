
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
