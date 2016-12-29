#' Converts TSS data into bed format
#' @param experimentName an object of class tssObject with tss data loaded
#' @param fileName the name of the BED file to be written
#' @export

setGeneric(
           name="tssToBed",
           def=function(experimentName, fileName="myTSS") {
               standardGeneric("tssToBed")
           }
    )

setMethod("tssToBed",
          signature(experimentName="tssObject", fileName="character"),

          function(experimentName, fileName="myTSS") {
              if (length(experimentName@tssData) == 0) {
                  stop("Slot @tssData is empty.\n\n Please load alignment files to your tssObject object.\n\n")
              }

              else {
                 experimentName@tssData -> my.TSS
                 cat("Converting", length(my.TSS), "GRanges objects to BED format...\n\n")

                  for (i in 1:length(my.TSS)) {
                      cat("Converting object number", i, "...\n\n")
                      my.TSS[[i]] -> this.GR
                      as(this.GR, "data.frame") -> this.df
                      nrow(this.df) -> my.len
                      rep(".", my.len) -> my_filler
                      cbind(as.character(this.df$seqnames), this.df$start, this.df$end, my_filler, my_filler, as.character(this.df$strand)) -> this.bed
                      c("chr", "start", "end", "ID", "score", "strand") -> colnames(this.bed)
                      paste(fileName, i, sep="_") -> this.name
                      paste(this.name, "bed", sep=".") -> this.name2
                      write.table(this.bed, file=this.name2, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
                  }

                  cat("\n")
                  cat(length(my.TSS), "BED files were written to your working directory.")
             }
          }
          )
