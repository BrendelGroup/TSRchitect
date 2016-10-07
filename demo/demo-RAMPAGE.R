### note: before you begin this demo, please copy the contents of /DATA/TSRchitect-demo/RAMPAGE into extdata/RAMPAGE.
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################
# Installation of TSRchitect from source
devtools::install_github("brendelgroup/TSRchitect")
# Loading the TSRchitect library
library(TSRchitect)
####################################################################################################
# initializing the tssExp object                                                                   
initializeExp("Example RAMPAGE experiment (Rice)", "riceRAMPAGE", "bam_filtered", isPairedEnd=TRUE)#
####################################################################################################
# loading the BAM files from the TSS profiling experiment (in this case, RAMPAGE)                     
importBam(riceRAMAPGE)                                                                             #
####################################################################################################
# converting BAM data into TSS information and attaching it to your tssExp object                  
bamToTSS(riceRAMPAGE)                                                                              #
####################################################################################################
# finding TSRs for a given dataset (in this case, the first one)
tsrFind(expName=riceRAMPAGE, tssNum=1, clustDist=20, nTSSs=3)				           #
####################################################################################################
