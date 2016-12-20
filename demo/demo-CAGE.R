### note: before you begin this demo, please copy the contents of /DATA/TSRchitect-demo/RAMPAGE into extdata/RAMPAGE.
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################

# Installation of TSRchitect from source
#devtools::install_github("brendelgroup/TSRchitect")
# Loading the TSRchitect library
library(TSRchitect)

# initializing the tss object:
initializeExp("Example human CAGE experiment", "CAGEhuman", c("extdata/CAGE"), isPairedEnd=FALSE) #

# setting the sample IDs and names for the datasets within the CAGE library
setSampleID(CAGEhuman, c("Example_1"), 1)

# loading the BAM files from the TSS profiling experiment (in this case, CAGE):
importBam(CAGEhuman)

# converting BAM data into TSS information and attaching it to your tssExp object:
bamToTSS(CAGEhuman)

# constructing a tss abundance matrix for the given dataset:
processTSS(expName=CAGEhuman, tssSet=1, writeTable=FALSE)

# finding TSRs for the given dataset:
tsrFind(expName=CAGEhuman, tssSet=1, nTSSs=5, clustDist=20, setToCluster="replicates", writeTable=FALSE)

####################################################################################################
#save(riceRAMPAGE, file="demo-CAGE-test.RData") #uncomment if you'd like to save a binary of the tssExp object to your working directory
