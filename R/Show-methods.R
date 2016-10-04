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
                  message("\nThe following bam files are loaded:\n")
                  for (i in 1:length(object@fileNames)) {
                      cat(paste(object@fileNames[i]),sep="\n")
                  }
              }
              else {
                  message("\nNo bam files loaded.\n")
              }
              if (length(object@tssData) > 0) {
                message("\nThere are ")
                cat(length(object@tssData))
                message(" tss datasets loaded in this object\n\n")
            }
              else {
                message("\nNo tss data loaded.\n")
            }
          }
          )


                  
