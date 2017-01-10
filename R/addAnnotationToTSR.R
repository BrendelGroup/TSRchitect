#' @title addAnnotationToTSR
#'
#' @description addAnnotationToTSR associates an identified promoter with a given gene, if found upstream and on the same strand within a selected distance.
#'
#' @param experimentName S4 object of class tssObject with filled data slots tsrData (and/or tsrDataMerged)
#' @param tsrSetType Specifies the type of TSR set to be processed.  Options are "replicates" or "merged".
#' @param tsrSet Number of the data set of type tsrSetType to be processed
#' @param upstreamDist - the maximum distance (in bp) upstream of the selected interval necessary to associate a TSR with a given annotation.
#' @param downstreamDist - the maximum distance (in bp) downstream of the selected interval (beginning with the CDS) to associate a TSR with a given annotation.
#' @param featureType - Does the annotation file have more than just genes in its 'type' column? Defaults to TRUE.
#' @param featureColumn - What is the name of the column in the annotation file containing the featureIDs? Defaults to 'ID'.
#' @param writeTable Specifies whether the output should be written to a tab-delimited file. Defaults to TRUE.
#'
#' @return addAnnotationToTSR adds featureID information to the tsrData(Merged) data frame and attaches it to its tssObject slot.
#'
#' @importFrom BiocGenerics start end
#' @importFrom GenomicRanges GRanges findOverlaps promoters
#' @importFrom IRanges IRanges
#' @importFrom utils write.table
#' @export


setGeneric(
           name="addAnnotationToTSR",
           def=function(experimentName, tsrSetType, tsrSet=1, upstreamDist, downstreamDist, featureType, featureColumn, writeTable=TRUE) {
               standardGeneric("addAnnotationToTSR")
    }
    )

setMethod("addAnnotationToTSR",
          signature(experimentName="tssObject", "character", "numeric", upstreamDist="numeric", downstreamDist="numeric", featureType="logical", featureColumn="character", "logical"),
          function(experimentName, tsrSetType, tsrSet, upstreamDist=1000, downstreamDist=200, featureType=TRUE, featureColumn="ID", writeTable=TRUE) {

             object.name <- deparse(substitute(experimentName))
             message("... addAnnotationToTSR ...")
             if (tsrSetType=="replicates") {
                 if (tsrSet>length(experimentName@tsrData)) {
                     stop("The value selected for tsrSet exceeds the number of slots in tsrData.")
                 }
                 outfname <- paste("TSRset-", tsrSet, sep="")
                 outfname <- paste(outfname, "txt", sep=".")
                 message("\nThe TSR set for TSS dataset ", tsrSet, " has been written to file ", outfname, "\nin your working directory.")
                 tsr.df <- experimentName@tsrData[[tsrSet]]
             }
             else if (tsrSetType=="merged") {
                 if (length(experimentName@tsrDataMerged)<1) {
                     stop("The @tsrDataMerged slot is currently empty. Please complete the merger before continuing.")
                 }
                 if (tsrSet>length(experimentName@tsrDataMerged)) {
                     stop("The value selected for tsrSet exceeds the number of slots in tsrDataMerged.")
                 }
                 if (tsrSet<length(experimentName@tssCountDataMerged)) {
                     outfname <- paste("TSRsetMerged-", tsrSet, sep="")
                     outfname <- paste(outfname, "txt", sep=".")
                     message("\nThe merged TSR set for TSS dataset ", tsrSet, " has been written to file ", outfname, "\nin your working directory.")
                 }
                 else { # "combined" case
                     outfname <- "TSRsetCombined.txt"
                     message("\nThe combined TSR set derived from all samples has been written to file ", outfname, "\nin your working directory.")
                 }
                 tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
             }
             else {
                 stop("Error: argument tsrSetType to addAnnotationToTSR() should be either \"replicates\" or \"merged\".")
             }

# ... loading the annotion GRange object and pulling out what we want:
             my.annot <- experimentName@geneAnnot
             if (length(my.annot)<1) {
                 stop("No annotation has been loaded to the tssObject. \nPlease run importAnnotation prior to using addAnnotationToTSR.")
             }
             if (featureType==TRUE) {
                 annot.gr <- my.annot[my.annot$type=="gene", ]
#                annot.gr <- my.annot[my.annot$type=="mRNA", ]
#VB COMMENT: ... this remains an open problem.  I think we want to specify two variable to be passed
#            into this function: 1) the "feature" tag we wish to pull out of my.annot, and 2) the column identifier in annot.gr to indicate what we'll use for annotation
             }
             else {
                 annot.gr <- my.annot
             }

             # ... creating a GRanges object from the data frame tsr.df:
             tsr.gr <- GRanges(seqnames=tsr.df$seq,
                               ranges = IRanges(start=tsr.df$start, end=tsr.df$end),
                               strand=tsr.df$strand)

# ... defining the regions of interest for annotation. Typically this would be the predicted promoter regions based on the gene annotion.
# promoters() takes the 5'-start at x of annot.gr, then creates an interval x-upstreamDist,x+downstreamDist-1 on the plus strand
# or x-downstreamDist+1,x+upstreamDist on the minus strand; in short - a potential promoter region.
             annot.extend <- promoters(annot.gr, upstream=upstreamDist, downstream=downstreamDist) #extending the annotations interval.

             idvec <- sprintf("annot.extend$%s",featureColumn)
             ID.vec <- eval(parse(text=idvec))
             my.OL <- findOverlaps(tsr.gr, annot.extend)
# ... my.OL is a hit list that indicates the overlaps between tsr.gr entries and annot.extend entries:
             OL.df <- as.data.frame(my.OL)

#... adding featureID to tsr.df:
#
             tsr.df$featureID <- NA #seeding the data frame (with NAs, which will represent no overlap after the next line of code is complete

	     # now replacing the NAs in the queryHits locations with the featureColumn entries of corresponding subject = annot.extend entries:

             tsr.df$featureID[OL.df$queryHits] <- ID.vec[OL.df$subjectHits]
             rownames(tsr.df) <- paste(tsr.df$seq, tsr.df$start, tsr.df$end, tsr.df$strand, sep=".") #adding the promoterIDs to the rows of the tsr data frame
	     # not sure whether we want a column or rownames with this ID ...

             if (writeTable=="TRUE") {
                 write.table(tsr.df, file=outfname, col.names=NA, row.names=TRUE, sep="\t", quote=FALSE)
                 message("\nThe updated TSR data have been written to file ", outfname, "\nin your working directory.")
             }
             #Update the record:
             if (tsrSetType=="replicates") {
                 tsr.df -> experimentName@tsrData[[tsrSet]]
             }
             else {
                 tsr.df -> experimentName@tsrDataMerged[[tsrSet]]
             }
             
             cat("Done. GeneIDs have been associated with adjacent TSRs and the data frame has been re-assigned to its slot.\n")
             cat("--------------------------------------------------------------------------------\n")
             assign(object.name, experimentName, envir = parent.frame())
             message(" Done.\n")
          }
          )
