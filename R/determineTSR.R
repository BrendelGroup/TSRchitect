#' determineTSR
#' Finds TSRs from a given sequence
#'
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param setToCluster - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#'
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="determineTSR",
           def=function(experimentName, setToCluster, tagCountThreshold, clustDist, writeTable=FALSE) {
               standardGeneric("determineTSR")
    }
    )

setMethod("determineTSR",
          signature(experimentName="tssObject", "character", "numeric", "numeric", "logical"),

          function(experimentName, setToCluster, tagCountThreshold=1, clustDist, writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... determineTSR ...")
             if (setToCluster=="replicates") {
                 iend <- length(experimentName@tssCountData)
                 experimentName@tsrData <- foreach(i=1:iend,.packages="TSRchitect") %dopar% findTSR(experimentName = experimentName, setToCluster="replicates", tssSet = i, tagCountThreshold=10, clustDist=20, writeTable = TRUE)
             }
	     else if (setToCluster=="merged") {
                 iend <- length(experimentName@tssCountDataMerged)
                 experimentName@tsrDataMerged <- foreach(i=1:iend,.packages="TSRchitect") %dopar% findTSR(experimentName = experimentName, setToCluster="merged", tssSet = i, tagCountThreshold=10, clustDist=20, writeTable = TRUE)

#Now we determine the TSS tag counts within the combined TSR set for each of the samples ...
                 experimentName@tsrDataMerged[[iend]] -> combinedTSRset

# ... going over each sample (index j):
                 for (j in 1:length(experimentName@tssCountData)) {
                    experimentName@tssCountData[[j]] -> this.tssSet
#... we are discarding counts below the tag count threshold tagCountThreshold:
                    this.tssSet <-this.tssSet[this.tssSet$nTSSs >= tagCountThreshold, ]
                    countv <- numeric(nrow(combinedTSRset))
                    if (nrow(this.tssSet) > 0) { # ... the following only makes sense if this.tssSet is non-empty

#... now considering each combined TSR in turn (index k):
                        lasttsrseq <- "null"
                        for (k in 1:nrow(combinedTSRset)) {
                            combinedTSRset[k,] -> this.tsr

                            if (this.tsr$seq != lasttsrseq) {
#... nothing to be counted if the current TSS seq does not match the current TSR seq:
                                that.tssSet <- this.tssSet[this.tssSet$seq == this.tsr$seq, ]
                                1 -> lbeg
                                0 -> count
                                lasttsrseq <- this.tsr$seq
                            }
                            if (nrow(that.tssSet) == 0) { # ... the following only makes sense if that.tssSet is non-empty
                                next
                            }
#... going over the TSS tags (index l) from sample j in order (first plus, then minus strand, position increasing)
                            for (l in lbeg:nrow(that.tssSet)) {
                               that.tssSet[l,] -> this.tss
# ... otherwise, we see how a given sample TSS position fits into the combined TSRs:
                               if (this.tsr$strand == this.tss$strand) {
                                  if (this.tsr$start  <= this.tss$TSS   &&   this.tss$TSS <= this.tsr$end) { 
                                     count + this.tss$nTSSs -> count
                                     if (l == nrow(that.tssSet)) {
                                         countv[k] <- count
                                     }
                                  }
                                  else if (this.tss$TSS  > this.tsr$start) {
                                     countv[k] <- count
                                     l -> lbeg
                                     0 -> count
                                     break
                                  }
                               }
                               else { # ... strand mismatch: go on to next comparison 
                                   if (this.tsr$strand == "+") { # ... need to go to next TSR
                                       countv[k] <- count
                                       l -> lbeg
                                       0  -> count
                                       break
                                   }
                               }
                            } # end for (l ...
                        } # end for (k ...
                    }
                    combinedTSRset$cname <- countv
                    colnames(combinedTSRset)[which(names(combinedTSRset) == "cname")] <- experimentName@sampleNames[j]

                 } # end for (j ...
                 if (writeTable=="TRUE") {
                     df.name <- paste(object.name, "TSR-tagcounts", sep="-")
                     df.name <- paste(df.name, "txt", sep=".")
                     write.table(combinedTSRset, file=df.name, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
                     message("\nThe TSR data have been written to file ", df.name, "\nin your working directory.")
                 }
             }
             else {
                 stop("Error: argument setToCluster to determineTSR() should be either \"replicates\" or \"merged\".")
             }
             cat("--------------------------------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
