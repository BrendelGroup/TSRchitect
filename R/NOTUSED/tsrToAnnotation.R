#' @title tsrToAnnotation
#'
#' @description tsrToAnnotation associates an identified promoter with a given gene, if found upstream and on the same strand within a selected distance.
#'
#' @param experimentName - an S4 object of class tssObject that contains information about the experiment.
#'
#' @param upstreamDist - the maximum distance (in bp) upstream of the selected interval necessary to associate a TSR with a given annotation.
#'
#' @param downstreamDist - the maximum distance (in bp) downstream of the selected interval (beginning with the CDS) to associate a TSR with a given annotation.
#'
#' @param featureType - Does the annotation file have more than just genes in its 'type' column? Defaults to TRUE.
#'
#' @param featureColumn - What is the name of the column in the annotation file containing the featureIDs? Defaults to 'ID'.
#'
#' @return tsrToAnnotation adds featureID information to the tsrData data frame and attaches it to the tssObject
#' experimentName
#'
#' @importFrom BiocGenerics start end
#' @importFrom GenomicRanges GRanges findOverlaps promoters
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom IRanges IRanges
#' @export

setGeneric(
    name="tsrToAnnotation",
    def=function(experimentName, upstreamDist, downstreamDist, featureType, featureColumn) {
        standardGeneric("tsrToAnnotation")
    }
    )

setMethod("tsrToAnnotation",
          signature(experimentName="tssObject", upstreamDist="numeric", downstreamDist="numeric", featureType="logical", featureColumn="character"),
          function(experimentName, upstreamDist=1000, downstreamDist=200, featureType=TRUE, featureColumn="ID") {
              experimentName.seq <- deparse(substitute(experimentName))
              message("... tsrToAnnotation ...")
              my.annot <- experimentName@geneAnnot
              if (length(my.annot)<1) {
                  stop("No annotation has been loaded to the tssObject. \nPlease run importAnnotation prior to using tsrToAnnotation.")
              }
              if (featureType==TRUE) {
                  annot.df <- my.annot[my.annot$type=="gene", ]
              }
              else {
                  annot.df <- my.annot
              }
              n.sets  <- length(experimentName@tsrData)
              tsr.set <- experimentName@tsrData[[n.sets]] #creating a GRanges object from the data frame in tsr.set. There's no straightforward way to do this without doing what follows, because the start and end columns need to be handled separately.
              df.plus <- tsr.set[tsr.set$strand=="+",]
              df.minus <- tsr.set[tsr.set$strand=="-",]
              gr1 <- GRanges(seqnames=df.plus$seq,  
                                      ranges = IRanges(
                                          start=df.plus$start,
                                          end=df.plus$end
                                          ),
                                      strand=df.plus$strand
                                      )
              gr2 <- GRanges(seqnames=df.minus$seq,
                                       ranges = IRanges(
                                           start=df.minus$start,
                                           end=df.minus$end
                                           ),
                                       strand=df.minus$strand
                                       )
              gr.combined <- c(gr1,gr2)
              gr.combined <- sortSeqlevels(gr.combined)
              tsr.gr <- sort(gr.combined) #The TSR data is now a GRanges object
              #extending the gene annotation upstream as specified by upstreamDist
              annot.extend <- promoters(annot.df, upstream=upstreamDist, downstream=downstreamDist) #extending the annotations interval.
              ID.vec <- annot.extend$featureColumn
              my.OL <- findOverlaps(tsr.gr, annot.extend) #finding which intervals overlap between the tsr data and the annotation (which has been extended upstream as specified.
              OL.df <- as.data.frame(my.OL)
              tsr.set$featureID <- NA #seeding the data frame (with NAs, which will represent no overlap after the next line of code is complete
              tsr.set$featureID[OL.df$queryHits] <- ID.vec[OL.df$subjectHits]
              rownames(tsr.set) <- paste(tsr.set$seq, tsr.set$start, tsr.set$end, tsr.set$strand, sep=".") #adding the promoterIDs to the rows of the tsr data frame
              n.merged <- length(experimentName@tsrDataMerged)
              n.slot <- n.merged + 1
              experimentName@tsrDataMerged[[n.slot]] <- tsr.set #(for now). I still think we should create a new, dedicated slot for merged data.
              cat("Done. GeneIDs have been associated with adjacent TSRs and the data frame has been re-assigned to its slot.\n")
              cat("--------------------------------------------------------------------------------\n")
              assign(experimentName.seq, experimentName, parent.frame())
              message(" Done.\n")
          }
          )
