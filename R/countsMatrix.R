#' Creates a data frame with a list of TSS positions and the corresponding counts from each
#' @param experimentName an S4 object of class tssObject that contains information about the experiment (including tssData)
#' @return Assigns a data frame to the expData slot of your tssObject
#' @importFrom gtools mixedsort
#' @export

setGeneric(
    name="countsMatrix",
    def=function(experimentName) {
        standardGeneric("countsMatrix")
    }
    )

setMethod("countsMatrix",
          signature(experimentName="tssObject"),
          function(experimentName) {
              object.name <- deparse(substitute(experimentName))

              message("\nTSS conversion is underway.")

              if (length(experimentName@tssData) == 0) {
                  stop("\nAll tssData slots are empty.\nYou must load your TSS data before proceeding with this function.")
              }
              n.slots <- length(experimentName@tssData)
              tss.list <- vector(mode="list", length=n.slots)
              names.vec <- vector(mode="character")

              for (i in 1:n.slots) {
                  paste("TSSset", i, sep="_") -> this.name
                  processTSS(experimentName, i, FALSE)
                  experimentName@expData[[i]] -> my.exp
                  my.exp -> tss.list[[i]]
                  paste(as.character(my.exp$chr), my.exp$CTSS, sep="-") -> my.names
                  paste(my.names, my.exp$strand, sep="_") -> these.names
                  c(names.vec, these.names) -> names.vec
              }
                  ctss.names <- unique(names.vec)
                  my.len <- length(ctss.names)
                  names.sorted <- mixedsort(ctss.names)
                  last.matrix <- matrix(NA, nrow=my.len, ncol=n.slots)
                  rownames(last.matrix) <- names.sorted
              for (j in 1:n.slots) {
                  tss.list[[j]] -> this.exp
                  paste(as.character(this.exp$chr), this.exp$CTSS, sep="-") -> exp.names.i
                  paste(exp.names.i, this.exp$strand, sep="_") -> exp.names
                  for (l in 1:length(ctss.names)) {
                      names.sorted[l] -> one.name
                      which(exp.names==one.name) -> this.ind #getting the index in the table so I can retrieve the counts info
                      if (length(this.ind)==0) {
                          0 -> this.string
                      }
                      else {
                          this.exp[this.ind[1], 3] -> this.string
                      }
                      if (is.na(this.string)==TRUE) {
                          stop("Here's a problem")
                      }
                      this.string -> last.matrix[l,j]
                  }
              }
                  last.frame <- as.data.frame(last.matrix)
                  colnames(last.frame) <- experimentName@sampleNames
                  experimentName@countsData <- last.frame
                  cat("Done. TSS abundance data was successfully added to your tssObject\n", object.name,"\".\n")
                  assign(object.name, experimentName, envir = parent.frame())
                  message(" Done.\n")
          }
          )
