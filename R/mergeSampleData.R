#' @title \strong{mergeSampleData}
#' @description \code{mergeSampleData} combines samples from multiple TSS
#' experiments into a single \linkS4class{GRanges} object
#'
#' @param experimentName an S4 object of class \emph{tssObject} that contains
#' information about the experiment.
#' @param n.cores the number of cores to be used for this job.
#' ncores=1 means serial execution of function calls (numeric)
#' @param tagCountThreshold the number of TSSs required at a given position
#' for it to be considered in sample data merging. (numeric)
#' Note: Merged data sets can become very large when the tagCountThreshold is
#'       set low (leading to inclusion of a lot of "noise" and resulting in
#'       long execution times in the necessary TSS position ordering step).
#'
#' @return tssCountData datasets are merged (according to the
#' \emph{sampleIDs}) and put in the tssCountDataMerged slot in the returned
#' \emph{tssObject}.
#'
#' @importFrom GenomicRanges as.data.frame
#' @importFrom gtools mixedorder
#' @importFrom utils flush.console
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' tssObjectExample <- mergeSampleData(experimentName=tssObjectExample,
#'                                     n.cores=1, tagCountThreshold=1)
#'
#' @note An example similar to the one provided can be found in
#' the vignette (/inst/doc/TSRchitect.Rmd).
#' @export
#' @rdname mergeSampleData-methods


setGeneric("mergeSampleData",
    function(experimentName, n.cores, tagCountThreshold)
    standardGeneric("mergeSampleData")
)

#' @rdname mergeSampleData-methods

setMethod("mergeSampleData",
          signature(experimentName="tssObject", "numeric", "numeric"),
          function(experimentName, n.cores=1, tagCountThreshold=1) {

              message("... mergeSampleData ...")
              if (length(experimentName@tssCountData)==0) {
                  stop("\nThe slot @tssCountData is empty.",
                       "Please run processTSS before proceeding with",
                       "this command.\n")
              }

              if (length(experimentName@sampleNames) < 1) {
                  stop("\nThe slot @sampleNames on your tssObject",
                       "object is empty. Please add sampleNames to",
                       "the object.\n")
              }

              if (length(experimentName@replicateIDs) < 1) {
                  stop("\nThe slot @replicateIDs on your tssObject",
                       "object is empty.\n",
                       "Please add replicateIDs to the object.\n")
              }

# mergeTSSdf() takes a combined data frame of TSS count data and collapses
# multiple entries for the same TSS position (from the different samples
# being merged) to a single entry with cumulated tag counts:
#
              mergeTSSdf <- function(cTSSdf) {
                mTSSdf <- data.frame(seq=character(nrow(cTSSdf)),
                                     TSS=numeric(nrow(cTSSdf)),
                                     strand=character(nrow(cTSSdf)),
                                     nTAGs=integer(nrow(cTSSdf)),
                                     stringsAsFactors=FALSE         )
                if (nrow(cTSSdf) <= 1) {
                  return(cTSSdf)
                }
                j = 1
                crow <- cTSSdf[1,]
                for (k in 2:nrow(cTSSdf)) {
                   nrow <- cTSSdf[k,]
                   if (nrow$TSS == crow$TSS  &  nrow$seq == crow$seq  &
                       nrow$strand == crow$strand) {
                     crow$nTAGs <- crow$nTAGs + nrow$nTAGs
                   } else {
                     mTSSdf[j,] <- crow
                     crow <- nrow
                     j <- j + 1
                   }
                   if ( k %% 25000 == 0) {
                     message("... done merging ", k, " of ",nrow(cTSSdf),
			     " TSS positions at\t", date())
		     flush.console()
		   }
                }
                mTSSdf[j,] <- crow
                return(mTSSdf[1:j,])
              }

              if (n.cores > 1) {
                BiocParallel::register(MulticoreParam(workers=n.cores),
                                                      default=TRUE)
              }

              rep.ids <- experimentName@replicateIDs
              uni.ids <- unique(rep.ids)
# ... we'll ignores samples with replicateID equal to zero as well as TSS
#     with less than tagCountThreshold support:
#
              uni.ids <- uni.ids[uni.ids > 0]
              uni.slots <- length(uni.ids)
              exp.data <- lapply(experimentName@tssCountData,
                                 function(df) {
                                   df[df$nTAGs >= tagCountThreshold,] })
              exp.list <- vector(mode="list")
              df.ind <- lapply(seq_along(uni.ids),
                               function(i) { which(rep.ids==i) })

              if (n.cores > 1) {
                exp.list <- bplapply(seq_along(df.ind), function(i) { 
                  this.df <- do.call(rbind, exp.data[df.ind[[i]]])
                  this.df <- this.df[order(this.df$seq, this.df$TSS),]
                  this.df <- this.df[mixedorder(this.df$seq),]
                  mergeTSSdf(this.df)                                  })
              } else {
                exp.list <- lapply(seq_along(df.ind), function(i) { 
                  this.df <- do.call(rbind, exp.data[df.ind[[i]]])
                  this.df <- this.df[order(this.df$seq, this.df$TSS),]
                  this.df <- this.df[mixedorder(this.df$seq),]
                  mergeTSSdf(this.df)                                  })
              }

# The following few lines merge the merged tssCountData into the last
# experimentName@tssCountDataMerged slot, representing the entire
# collection of TSS tag counts in the experiment:
#
              n.slots <- length(uni.ids) + 1
              my.df <- do.call(rbind, exp.list)
              my.df <- my.df[order(my.df$seq, my.df$TSS),]
              my.df <- my.df[mixedorder(my.df$seq),]
              exp.list[[n.slots]] <- mergeTSSdf(my.df)
              experimentName@tssCountDataMerged <- exp.list

              message("\n... the TSS expression data have been merged",
                    "\nand added to the tssObject object.\n")
              message("------------------------------------------------------\n")
              message(" Done.\n")

              return(experimentName)
          }
          )
