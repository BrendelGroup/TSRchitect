# Show methods for TSRchitect
      
setMethod("show",
          signature(object="tssExp"),
          function(object) {
              cat("############################################\n")
              cat("############## TSRchitect ##################\n")
              cat("       Object of Class 'tssExp'\n")
              cat("############################################\n\n")
              cat("Title of experiment:", object@title,"\n\n")
              cat("The TSS data were specified to be \"", object@dataType, "\"\n")
              if (length(object@fileNames) > 0) {
                  cat("and to be contained in the following files:\n")
                  for (i in 1:length(object@fileNames)) {
                      cat(paste(object@fileNames[i]),sep="\n")
                  }
              }
              else {
                  stop("\nNo *.bam files were found.  Please check.\n")
              }
              if (length(object@tssData) > 0) {
                  expName <- deparse(substitute(object))
                  cat("\n")
                  cat(length(object@tssData))
                  cat(" TSS datasets were loaded into the tssExp object \"", expName, "\".\n")
              }
              else {
                  cat("\nNo tss data have been loaded.\n")
              }
          }
          )
