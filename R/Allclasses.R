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

tssObject <- function(title=NA, dataType=NA, bamData=NA) {

    new.tssObj <- new("tssObject",
                      title=title,
                      dataType=dataType,
                      bamData=bamData)

    if (!(is.na(title))) {
        new.tssObj@title <- title
    }

    if (!(is.na(dataType))) {
        new.tssObj@dataType <- dataType
    }

    if (!(is.na(bamData))) {
        new.tssObj@bamData <- bamData
    }

    message("\nA new tssObject has been created in your workspace.\n")
    return(new.tssObj)
}
                      
