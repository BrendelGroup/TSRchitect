#' tsrFind
#' Finds TSRs from a given sequence
#' @param experimentName - a S4 object of class tssObject containing information in slot tssTagData
#' @param tssSet - number of the dataset to be analyzed
#' @param tagCountThreshold - number of TSSs required at a given position
#' @param clustDist - maximum distance of TSSs between two TSRs (in base pairs)
#' @param setToCluster - specifies the set to be clustered. Options are "replicates" or "merged".
#' @param writeTable - specifies whether the output should be written to a table. (logical)
#' @return creates a list of GenomicRanges containing TSR positions in slot 'tsrData' on your tssObject object
#' @export

setGeneric(
           name="tsrFind",
           def=function(experimentName, tssSet=1, tagCountThreshold, clustDist, setToCluster, writeTable=FALSE) {
               standardGeneric("tsrFind")
    }
    )

setMethod("tsrFind",
          signature(experimentName="tssObject", "numeric", "numeric", "numeric", "character", "logical"),

          function(experimentName, tssSet, tagCountThreshold=1, clustDist, setToCluster, writeTable=FALSE) {
             object.name <- deparse(substitute(experimentName))

             message("... tsrFind ...")
             if (setToCluster=="replicates") {
                 if (tssSet>length(experimentName@tssCountData)) {
                     stop("The value selected for tssSet exceeds the number of slots in tssTagData.")
                 }

                 tss.mat <- experimentName@tssCountData[[tssSet]]
                 tsr.list <- tsrCluster(tss.mat, minNbrTSSs=tagCountThreshold, minDist=clustDist)
                 tsr.DF <- tsrToDF(tsr.list)

                 if (writeTable=="TRUE") {
                     df.name <- paste("TSRset-", tssSet, sep="")
                     df.name <- paste(df.name, "txt", sep=".")
                     write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                     message("\nThe TSR set for TSS dataset ", tssSet, " has been written to file ", df.name, "\nin your working directory.")
                 }

                 experimentName@tsrData[[tssSet]] <- tsr.DF
                 cat("\n... the TSR data frame for dataset ", tssSet, " has been successfully added to\ntssObject object \"", object.name, "\"\n")
              }

              else if (setToCluster=="merged") {
                  if (length(experimentName@tssCountDataMerged)<1) {
                      stop("The @tssCountDataMerged slot is currently empty. Please complete the merger before continuing.")
                  }

                  tsr.list <- vector(mode="list")
                  for (i in 1:length(experimentName@tssCountDataMerged)) {
                      tss.mat <- experimentName@tssCountDataMerged[[i]]
                      my.tsr <- tsrCluster(tss.mat, minNbrTSSs=tagCountThreshold, minDist=clustDist)
                      tsr.DF <- tsrToDF(my.tsr)
                      tsr.list[[i]] <- tsr.DF

                      if (writeTable=="TRUE") {
                          if (i < length(experimentName@tssCountDataMerged)) {
                              df.name <- paste("TSRsetMerged-", i, sep="")
                              df.name <- paste(df.name, "txt", sep=".")
                              message("\nThe merged TSR set for TSS dataset ", i, " has been written to file ", df.name, "\nin your working directory.")
                          }
                          else {
                              df.name <- "TSRsetCombined.txt"
                              message("\nThe combined TSR set derived from all samples has been written to file ", df.name, "\nin your working directory.")
                          }
                          write.table(tsr.DF, file=df.name, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
                      }
                  }

                  experimentName@tsrDataMerged <- tsr.list
                  cat("\n... merged TSR data frames have been successfully added to\ntssObject object \"", object.name, "\"\n")

#Now we determined the TSS tag counts within the combined TSR set for each of the samples ...
                  length(experimentName@tssCountDataMerged) -> i
                  experimentName@tsrDataMerged[[i]] -> combinedTSRset

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
                  stop("Error: argument setToCluster to tsrFind() should be either \"replicates\" or \"merged\".")
              }
              cat("--------------------------------------------------------------------------------\n")
              assign(object.name, experimentName, envir = parent.frame())
              message(" Done.\n")
          }
          )
