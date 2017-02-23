#' @title \strong{addAnnotationToTSR}
#' @description \code{addAnnotationToTSR} associates an identified promoter
#' with a given gene, if found upstream and on the same strand within
#' a specified range.
#'
#' @param experimentName an object of class \emph{tssObject} with occupied
#' data slots \emph{@@tsrData} (and/or \emph{@@tsrDataMerged}).
#' The \emph{tssObject} must alrady have an annotation attached to the slot
#' \emph{@@annotation}, which is provided by either
#' \code{\link{importAnnotationExternal}} or
#' \code{\link{importAnnotationHub}}.
#' @param tsrSetType Specifies the type of TSR set to be processed.
#' Options are "replicates" or "merged".
#' @param tsrSet Number of the data set (of type \emph{tsrSetType}) to be
#' processed. (numeric)
#' @param upstreamDist the maximum distance (in bp) upstream of the selected
#' interval necessary to associate a TSR with a given annotation. (numeric)
#' @param downstreamDist the maximum distance (in bp) downstream of the start
#' of the selected interval to associate a TSR with a given annotation.
#' (numeric)
#' @param feature Specifies the feature to be used for annotation
#' (typically "gene" [default] or "mRNA" for GFF3 input); set to "all" if all
#' annotations from the input are to be used. (character)
#' @param featureColumnID Name of the column identifier in the
#' \linkS4class{GRanges} annotation object. This should be "ID" (default) for
#' GFF3 input or "name" for bed input. (character)
#' @param writeTable logical, specifying whether the output should be written
#' to a tab-delimited file. Defaults to TRUE.
#`
#' @return addAnnotationToTSR adds feature annotation to the (merged)
#' \emph{@@tsrData} data frame and returns the updated \emph{tssObject}.
#'
#' @importFrom BiocGenerics start end
#' @importFrom GenomicRanges GRanges findOverlaps promoters
#' @importFrom IRanges IRanges
#' @importFrom utils write.table
#'
#' @examples
#' load(system.file("extdata", "tssObjectExample.RData",
#' package="TSRchitect"))
#' addAnnotationToTSR(experimentName=tssObjectExample, tsrSetType="merged",
#' tsrSet=1, upstreamDist=1000, downstreamDist=200, feature="transcript",
#' featureColumnID="ID", writeTable=FALSE)
#' #if the object attached to @@annotation is a gff/gff3 file
#'
#' @note An example similar to the this one can be found
#' in the vignette (/inst/doc/TSRchitect.Rmd)
#'
#' @export


setGeneric("addAnnotationToTSR",
    function(experimentName, tsrSetType, tsrSet=1, upstreamDist,
             downstreamDist, feature, featureColumnID, writeTable=TRUE)
    standardGeneric("addAnnotationToTSR")
)

setMethod("addAnnotationToTSR",
          signature(experimentName="tssObject", "character", "numeric",
                    upstreamDist="numeric", downstreamDist="numeric",
                    feature="character", featureColumnID="character",
                    "logical"),
          function(experimentName, tsrSetType, tsrSet, upstreamDist=1000,
                   downstreamDist=200, feature="gene", featureColumnID="ID",
                   writeTable=TRUE) {

              object.name <- deparse(substitute(experimentName))
              message("... addAnnotationToTSR ...")
              if (tsrSetType=="replicates") {
                  if (tsrSet>length(experimentName@tsrData)) {
                      stop("The value selected for tsrSet exceeds the number ",
                           "of slots in tsrData.")
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
                      stop("The @tsrDataMerged slot is currently empty.\n",
                          "Please complete the merger before continuing.")
                  }
                  if (tsrSet>length(experimentName@tsrDataMerged)) {
                      stop("The value selected for tsrSet exceeds the ",
                           "number of slots in tsrDataMerged.")
                  }
                  if (tsrSet<length(experimentName@tssCountDataMerged)) {
                      outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                      outfname <- paste(outfname, "txt", sep=".")
                      message("\nThe merged TSR set for TSS dataset ",
                              tsrSet, " has been written to file ", outfname,
                              "\nin your working directory.")
                  }
                  else { # "combined" case
                      outfname <- "TSRsetCombined.txt"
                      message("\nThe combined TSR set derived from all samples",
                              " has been written to file ", outfname,
                              "\nin your working directory.")
                  }
                  tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
              }
              else {
                  stop("Error: argument tsrSetType to addAnnotationToTSR()",
                       "should be either \"replicates\" or \"merged\".")
              }

#  ... loading the annotion GRange object and pulling out what we want:
              allAnnotation <- experimentName@annotation
              if (length(allAnnotation)<1) {
                  stop("No annotation has been loaded to the tssObject.",
                       "\n\nPlease import an annotation prior to using ",
                       "addAnnotationToTSR.")
              }
#  ... pulling out the selected feature from the annotation GRanges
#  object (if the required, GFF3-standard, "type" column exists):
              if ( feature == "all"  || is.na(match("type",
                       names(allAnnotation@elementMetadata@listData))) ) {
                  annot.gr <- allAnnotation
              }
              else {
                  annot.gr <- allAnnotation[allAnnotation$type==feature, ]
              }

              # ... creating a GRanges object from the data frame tsr.df:
              tsr.gr <- GRanges(seqnames=tsr.df$seq,
                                ranges = IRanges(start=tsr.df$start,
                                    end=tsr.df$end),
                                strand=tsr.df$strand)

#  ... defining the regions of interest for annotation.
#  Typically this would be the predicted promoter regions based on
#  the gene annotion.
#  promoters() takes the 5'-start at x of annot.gr, then creates an interval
#  x-upstreamDist,x+downstreamDist-1 on the plus strand
#  or x-downstreamDist+1,x+upstreamDist on the minus strand; in short
#-  a potential promoter region.
              regionOfInterest <- promoters(annot.gr, upstream=upstreamDist,
                                            downstream=downstreamDist)

              idvec <- sprintf("regionOfInterest$%s",featureColumnID)
              ID.vec <- eval(parse(text=idvec))
              overlapHitList <- findOverlaps(tsr.gr, regionOfInterest)
#  ... overlapHitList is a hit list that indicates the overlaps
#  between tsr.gr entries and regionOfInterest entries:
              overlap.df <- as.data.frame(overlapHitList)

#...  adding featureID to tsr.df:
#
              tsr.df$featureID <- NA #seeding the data frame (with NAs,
#  which will represent no overlap after the next line of code is complete

#  now replacing the NAs in the queryHits locations with the featureColumnID
#  entries of corresponding subject = regionOfInterest entries:

              tsr.df$featureID[overlap.df$queryHits] <-
              ID.vec[overlap.df$subjectHits]
              rownames(tsr.df) <- paste(tsr.df$seq, tsr.df$start, tsr.df$end,
                                        tsr.df$strand, sep=".")
#adding  the promoterIDs to the rows of the tsr data frame

              if (writeTable=="TRUE") {
                  write.table(tsr.df, file=outfname, col.names=NA,
                              row.names=TRUE,  sep="\t", quote=FALSE)
                  message("\nThe updated TSR data have been written to ",
                          "file ", outfname, " in your working directory.")
              }
              #Update the record:
              if (tsrSetType=="replicates") {
                  tsr.df -> experimentName@tsrData[[tsrSet]]
              }
              else {
                  tsr.df -> experimentName@tsrDataMerged[[tsrSet]]
              }

              cat("Done. GeneIDs have been associated with adjacent TSRs.\n")
              cat("---------------------------------------------------------\n")
              message(" Done.\n")
              return( experimentName)
          }
          )
