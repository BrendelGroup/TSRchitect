### note: before you begin this demo, please copy the contents of /DATA/TSRchitect-demo/RAMPAGE into extdata/RAMPAGE.
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################

# Installation of TSRchitect from source
#devtools::install_github("brendelgroup/TSRchitect")
# Loading the TSRchitect library
library(TSRchitect)

# initializing the tssExp object:
initializeExp("Example RAMPAGE experiment (Rice)", "riceRAMPAGE", "bam_filtered", isPairedEnd=TRUE)

# setting the sample IDs and names for the datasets within the RAMPAGE library
setSampleID(riceRAMPAGE, c("DL1_1","DL1_2", "DL2_1", "DL2_2", "DL2_3"), c(1,1,2,2,2))

# loading the BAM files from the TSS profiling experiment (in this case, RAMPAGE):
importBam(riceRAMPAGE)

# converting BAM data into TSS information and attaching it to your tssExp object:
bamToTSS(riceRAMPAGE)

# constructing a tss abundance matrix for the given dataset:
processTSS(expName=riceRAMPAGE, tssSet=1, writeTable=FALSE)
processTSS(expName=riceRAMPAGE, tssSet=2, writeTable=FALSE)
processTSS(expName=riceRAMPAGE, tssSet=3, writeTable=FALSE)
processTSS(expName=riceRAMPAGE, tssSet=4, writeTable=FALSE)
processTSS(expName=riceRAMPAGE, tssSet=5, writeTable=FALSE)

# finding TSRs for the given dataset:
tsrFind(expName=riceRAMPAGE, tssSet=1, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)
tsrFind(expName=riceRAMPAGE, tssSet=2, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)
tsrFind(expName=riceRAMPAGE, tssSet=3, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)
tsrFind(expName=riceRAMPAGE, tssSet=4, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)
tsrFind(expName=riceRAMPAGE, tssSet=5, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)

# merging tssData objects according to the replicate info
mergeTSS(expName=riceRAMPAGE)

# merging exprData objects according to the replicate info
mergeExpr(expName=riceRAMPAGE)

# finding TSRs from the merged datasets
tsrFind(expName=riceRAMPAGE, tssSet=1, nTSSs=3, clustDist=20, setToCluster="merged", writeTable=FALSE)

####################################################################################################
save(riceRAMPAGE, file="demo-RAMPAGE-test.RData") #uncomment if you'd like to save a binary of the tssExp object to your working directory
