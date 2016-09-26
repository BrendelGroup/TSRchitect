# Show methods for TSRchitect
      
setMethod("show",
          signature(object= "tssExp"),
          function(object) {
              cat("############################################\n")
              cat("############## TSRchitect ##################\n")
              cat("       Object of Class 'tssExp'\n")
              cat("############################################\n")
              cat("Title of Experiment:", object@title,"\n")
              cat("TSS data is:", object@dataType, "\n")
              if (length(object@fileNames) > 0) {
                  message("The following bam files are loaded:\n")
                  for (i in 1:length(object@fileNames)) {
                      cat(paste(object@fileNames[i]))
                  }
              }
              else {
                  message("No bam files loaded.\n") }
              cat(length(object@bamData)," files are loaded in the @bamData slot.\n")
          }
          )


                  
