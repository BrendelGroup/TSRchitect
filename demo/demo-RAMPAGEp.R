### note: before you begin this demo, please copy the contents of /DATA/TSRchitect-demo/RAMPAGE into extdata/RAMPAGE.
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################

# Installation of TSRchitect from source
#devtools::install_github("brendelgroup/TSRchitect")
# Loading the TSRchitect library
library(TSRchitect)
library(doParallel)

# setting up a cluster of 6 nodes:
mycl <- makeCluster(6,type="FORK",outfile="")
registerDoParallel(mycl)

# initializing the tssExp object:
#initializeExp("Example RAMPAGE experiment (Rice)", "riceRAMPAGE", "bam_filtered", isPairedEnd=TRUE)
initializeExp("IRBB7 RAMPAGE experiment (Rice)", "riceRAMPAGE", "/scratch/OSseqPRJ/RAMPAGE/TSRchitect_test/IRBB7Aln", isPairedEnd=TRUE)

# setting the sample IDs and names for the datasets within the RAMPAGE library
setSampleID(riceRAMPAGE, c("DL1_1","DL1_2", "DL2_1", "DL2_2", "DL2_3"), c(1,1,2,2,2))

# loading the BAM files from the TSS profiling experiment (in this case, RAMPAGE):
importBam(riceRAMPAGE)

# converting BAM data into TSS information and attaching it to your tssExp object:
bamToTSS(riceRAMPAGE)

# constructing a tss abundance matrix for all datasets in parallel:
system.time(  riceRAMPAGE@expData <- foreach(i=1:5,.packages="TSRchitect") %dopar% processTSSp(experimentName = riceRAMPAGE, tssSet = i, writeTable = FALSE)  )


system.time(  mytestM@tsrData <- foreach(i=1:5,.packages="TSRchitect") %dopar% tsrFindP(experimentName = mytestM, tssSet = i, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable = TRUE)  )


# finding TSRs for the given dataset:
system.time(  riceRAMPAGE@tsrData <- foreach(i=1:5,.packages="TSRchitect") %dopar% tsrFindP(experimentName = riceRAMPAGE, tssSet = i, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable = FALSE)  )

# merging tssData objects according to the replicate info:
mergeTSS(experimentName=riceRAMPAGE)

# merging exprData objects according to the replicate info:
mergeExpr(experimentName=riceRAMPAGE)

# finding TSRs from the merged datasets:
tsrFind(experimentName=riceRAMPAGE, tssSet=1, nTSSs=3, clustDist=20, setToCluster="merged", writeTable=FALSE)

####################################################################################################
#save(riceRAMPAGE, file="demo-RAMPAGE-test.RData") #uncomment if you'd like to save a binary of the tssExp object to your working directory
