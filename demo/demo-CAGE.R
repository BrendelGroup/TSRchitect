### note: before you begin this test, please copy the contents of /DATA/TSRchitect-demo/ into extdata/CAGE
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################
library(TSRchitect)
# initializing the tssExp object                                                                   
initializeExp("Example human CAGE experiment", "CAGEhuman", c("extdata/CAGE"), isPairedEnd=FALSE) #
####################################################################################################
# loading the BAM files from the TSS profiling experiment (in this case, CAGE)                     
importBam(CAGEhuman)                                                                               #
####################################################################################################
# converting BAM data into TSS information and attaching it to your tssExp object                  
bamToTSS(CAGEhuman)                                                                                #
####################################################################################################
# finding TSRs for a given chromosome
tsrFind(expName=CAGEhuman, tssNum=1, clustDist=20, nTSSs=5)
####################################################################################################
