# Show methods for TSRchitect

setMethod("show",
          signature(object="tssObject"),
          function(object) {
              cat("############################################\n")
              cat("############## TSRchitect ##################\n")
              cat("     S4 Object of Class 'tssObject'\n")
              cat("############################################\n\n")
              if (is.na(object@title)==FALSE) {
                  cat("Title of experiment:", object@title,"\n\n")
              }
              else {
                  cat("No experiment title has been added.\n\n")
              }
              if (length(is.na(object@fileNamesBAM))>0) {
                  cat("There are ", length(object@fileNamesBAM))
                  cat(".bam datasets loaded, as follows:\n")
                  for (i in 1:length(object@fileNamesBAM)) {
                      cat(paste(object@fileNamesBAM[i]), sep="\n")
                  }
                  cat("\n")
              }
              else {
                  cat("\nNo *.bam datasets have been added. \n")
              }
              if (length(is.na(object@fileNamesBED))>0) {
                  cat("There are ", length(object@fileNamesBED))
                  cat(".bed datasets loaded, as follows:\n")
                  for (i in 1:length(object@fileNamesBED)) {
                      cat(paste(object@fileNamesBED[i]), sep="\n")
                  }
                  cat("\n")
              }
              else {
                  cat("\nNo *.bed datasets have been added. \n")
              }
              if (length(is.na(object@sampleNames))>0) {
                  cat("\nThe names of the datasets are:\n")
                  for (i in 1:length(object@sampleNames)) {
                      cat(paste(object@sampleNames[i]), sep="\n")
                  }
              }
              else {
                  cat("\nDataset names have not been provided.\n")
              }
              if (length(object@tssTagData) > 0) {
                  object.name <- as.character(deparse(substitute(object)))
                  cat("\n")
                  cat(length(object@tssTagData))
                  cat(" replicate TSS datasets were loaded into the",
                      "tssObject")
                  cat(".\n")
              }
              else {
                  cat("\nNo TSS data have been loaded.\n")
              }
              if (length(object@tssCountDataMerged) > 0) {
                  object.name <- as.character(deparse(substitute(object)))
                  cat("\n")
                  cat("Replicate datasets have been merged on the",
                      "tssObject")
                  cat(".\n")
              }
              else {
                  cat("\nTSS replicate datasets have not been merged.\n")
              }
             if (length(object@tsrData) > 0) {
                 object.name <- as.character(deparse(substitute(object)))
                 cat("\nTSRs have been identified from",
                     length(object@tsrData), "replicate datasets. \n")
              }
             else {
                 cat("\nTSRs have not been identified from replicate",
                    "datasets. \n")
             }
             if (length(object@tsrDataMerged) > 0) {
                 object.name <- as.character(deparse(substitute(object)))
                 cat("\nTSRs have been identified from",
                     length(object@tsrDataMerged), "merged datasets. \n")
              }
             else {
                 cat("\nTSRs have not been identified from",
                     "merged datasets.\n")
             }
          }
          )
