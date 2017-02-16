#' @title \strong{addTagCountsToTSR}
#' @description \code{addTagCountsToTSR} adds a matrix of tag counts to a set
#' of identified TSRs
#' 
#' @param experimentName a S4 object of class \emph{tssObject} containing
#' information in slot \emph{@tssTagData}
#' @param tsrSetType specifies the set to be written to file.
#' Options are "replicates" or "merged". (character)
#' @param tsrSet number of the dataset to be processed, where 1 corresponds
#' to the first slot, and so on. (numeric)
#' @param tagCountThreshold number of TSSs required at a given position
#' (numeric)
#' @param writeTable specifies whether the output should be written to
#' a table. (logical)
#' 
#' @return a matrix of tag counts (where the number of columns will equal
#' the number of replicates in the sample) is appended to the data frame
#' of the selected set of identified TSRs.
#' 
#' @importFrom utils write.table
#' 
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData", package="TSRchitect"))   
#' addTagCountsToTSR(experimentName=tssObjectExample, tsrSetType="merged",
#' tsrSet=1, tagCountThreshold=25, writeTable=FALSE)
#' 
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd)
#' 
#' @export

setGeneric(
           name="addTagCountsToTSR",
    def=function(experimentName, tsrSetType, tsrSet=1,
        tagCountThreshold=1, writeTable=TRUE) {
               standardGeneric("addTagCountsToTSR")
    }
    )

setMethod("addTagCountsToTSR",
          signature(experimentName="tssObject", "character", "numeric",
                    "numeric", "logical"),

          function(experimentName, tsrSetType, tsrSet, tagCountThreshold,
                   writeTable=TRUE) {
             object.name <- deparse(substitute(experimentName))

             message("... addTagCountsToTSR ...")
             if (tsrSetType=="replicates") {
                 if (tsrSet>length(experimentName@tsrData)) {
                     stop("The value selected for tsrSet exceeds the ",
                           "number of slots in tsrData.")
                 }
                 outfname <- paste("TSRset-", tsrSet, sep="")
                 outfname <- paste(outfname, "txt", sep=".")
                 message("\nThe TSR set for TSS dataset ", tsrSet,
                         " has been written to file ", outfname,
                         "\nin your working directory.")
                 tsr.df <- experimentName@tsrData[[tsrSet]]
             }
             else if (tsrSetType=="merged") {
                 if (length(experimentName@tsrDataMerged)<1) {
                     stop("The @tsrDataMerged slot is currently empty.",
                          "Please complete the merger before continuing.")
                 }
                 if (tsrSet>length(experimentName@tsrDataMerged)) {
                     stop("The value selected for tsrSet exceeds the ",
                          "number of slots in tsrDataMerged.")
                 }
                 if (tsrSet<length(experimentName@tssCountDataMerged)) {
                     outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                     outfname <- paste(outfname, "txt", sep=".")
                     message("\nThe merged TSR set for TSS dataset ", tsrSet,
                             " has been written to file ", outfname,
                             "\nin your working directory.")
                 }
                 else { # "combined" case
                     outfname <- "TSRsetCombined.txt"
                     message("\nThe combined TSR set derived from all samples ",
                             "has been written to file ", outfname,
                             "\nin your working directory.")
                 }
                 tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
             }
             else {
                 stop("Error: argument tsrSetType to addTagCountsToTSR() ",
                      "should be either \"replicates\" or \"merged\".")
             }

#Now we determine the TSS tag counts within the
#current TSR set for each of the samples ...
             tsr.df -> currentTSRset

# ... going over each sample (index j):
             for (j in 1:length(experimentName@tssCountData)) {
                if (experimentName@replicateIDs[j] == 0) {
 # ... samples with replicateID equal to zero are ignored
                    next
                }
                experimentName@tssCountData[[j]] -> this.tssSet
#... we are discarding counts below the tag count threshold tagCountThreshold:
                this.tssSet <-this.tssSet[this.tssSet$nTSSs >= tagCountThreshold, ]
                countv <- numeric(nrow(currentTSRset))
                if (nrow(this.tssSet) > 0) {
 # ... the following only makes sense if this.tssSet is non-empty

#... now considering each current TSR in turn (index k):
                    lasttsrseq <- "null"
                    for (k in 1:nrow(currentTSRset)) {
                        currentTSRset[k,] -> this.tsr

                        if (this.tsr$seq != lasttsrseq) {
#... nothing to be counted if the current TSS seq does not
# match the current TSR seq:
                            that.tssSet <- this.tssSet[this.tssSet$seq 
                                            == this.tsr$seq, ]
                            1 -> lbeg
                            0 -> count
                            lasttsrseq <- this.tsr$seq
                        }
                        if (nrow(that.tssSet) == 0) { 
# ... the following only makes sense if that.tssSet is non-empty
                            next
                        }
#... going over the TSS tags (index l) from sample j in order
# (first plus, then minus strand, position increasing)
                        for (l in lbeg:nrow(that.tssSet)) {
                           that.tssSet[l,] -> this.tss
# ... otherwise, we see how a given sample TSS position fits
# into the current TSRs:
                           if (this.tsr$strand == this.tss$strand) {
                              if (this.tsr$start  <= this.tss$TSS   &&  
                                  this.tss$TSS <= this.tsr$end) { 
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
                           else { # ... strand mismatch:
# go on to next comparison 
                               if (this.tsr$strand == "+") { 
# ... need to go to next TSR
                                   countv[k] <- count
                                   l -> lbeg
                                   0  -> count
                                   break
                               }
                           }
                        } # end for (l ...
                    } # end for (k ...
                }
                currentTSRset$cname <- countv
                colnames(currentTSRset)[which(names(currentTSRset) == 
                "cname")] <- experimentName@sampleNames[j]

             } # end for (j ...
             if (writeTable=="TRUE") {
                 write.table(currentTSRset, file=outfname, col.names=TRUE,
 row.names=FALSE, sep="\t", quote=FALSE)
                 message("\nThe updated TSR data have been written to file ",
 outfname, "\nin your working directory.")
             }
             #Update the record:
             if (tsrSetType=="replicates") {
                 currentTSRset -> experimentName@tsrData[[tsrSet]]
             }
             else {
                 currentTSRset -> experimentName@tsrDataMerged[[tsrSet]]
             }
             
             cat("---------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
