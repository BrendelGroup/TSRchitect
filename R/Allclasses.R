setClass(Class="tssExp",
         representation(
             title = "character",
             fileNames = "character",
             dataType = "character",
             bamData = "list",
             tssData = "GRangesList",
             expData = "data.frame",
             tsrData = "list"
             ),
         prototype(
             title =NA_character_,
             fileNames = NA_character_,
             dataType = NA_character_,
             bamData = list(),
             tssData = GRangesList(),
             expData = data.frame(),
             tsrData = list()
             ),
         )
            
