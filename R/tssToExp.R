#' Creates a data frame with a list of TSS positions and the corresponding counts from each
#' @param expName an S4 object of class tssExp that contains information about the experiment (including tssData)
#' @return Assigns a data frame to the expData slot of your tssExp object
#' @export 

setGeneric(
    name="tssToExp",
    def=function(expName) {
        standardGeneric("tssToExp")
    }
    )

setMethod("tssToExp",
          signature(expName="tssExp"),
          function(expName) {

              object.name <- deparse(substitute(expName))

              message("\nTSS conversion is underway.")

              if (length(expName@tssData) == 0) {

                  stop("\nAll tssData slots are empty.\nYou must load your TSS data before proceeding with this function.")
              }

              n.slots <- length(expName@tssData)

              tss.list <- vector(mode="list")

              for (i in 1:n.slots) {
                  paste("TSSset", i, sep="_") -> this.name
                  acquireTSS(expName, i) -> tss.set
                  expressionCTSS(tss.set, writeDF=FALSE) -> tss.mat
                  tss.mat -> tss.list[[this.name]]

              }
              my.exp <- tss.list[[1]]
              my.names <- paste(as.character(my.exp$chr), my.exp$CTSS, my.exp$strand, sep="-")

              for (j in 2:length(tss.list)) {
                  tss.list[[j]] -> this.exp
                  paste(as.character(this.exp$chr), this.exp$CTSS, my.exp$strand) -> these.names
                  c(my.names, these.names) -> my.names

              }
              ctss.names <- unique(my.names)
              my.len <- length(ctss.names)
              last.matrix <- matrix(NA, nrow=my.len, ncol=5)
              rownames(last.matrix) <- ctss.names

              for (k in 1:length(tss.list[[k]])) {
                  tss.list[[k]] -> this.exp
                  paste(as.character(this.exp$chr), this.exp$CTSS, my.exp$strand) -> exp.names
#                  data.frame(names=ctss.names,nTSSs=

                  for (l in 1:length(ctss.names)) {
                      ctss.names[l] -> one.name
                      which(exp.names==one.name) -> this.ind #getting the index in the table so I can retrieve the counts info
                      this.exp[this.ind, 2] -> this.string
                      as.character(this.string[1]) -> my.chr
                      as.numeric(as.character(this.string[2])) -> my.pos
                      as.numeric(as.character(this.string[3])) -> nTSS
                      as.numeric(as.character(this.string[4])) -> my.strand
                      c(my.chr, my.pos, nTSS, my.strand) -> last.matrix[l, ]
                      
              }

                  expName@expData <- tss.list
                  cat("\nTSS abundance data was successfully added to your tssExp object.\n")
                  return(tss.list)
 #            assign(object.name, expName, envir = parent.frame())              

          }
          )
