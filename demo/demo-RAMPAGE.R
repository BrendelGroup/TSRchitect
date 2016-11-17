### note: before you begin this demo, please copy the contents of /DATA/TSRchitect-demo/RAMPAGE into extdata/RAMPAGE.
### This demo should be run from the top-level directory of TSRchitect (TSRchitect/)
####################################################################################################
# Installation of TSRchitect from source
#devtools::install_github("brendelgroup/TSRchitect")
# Loading the TSRchitect library
library(TSRchitect)
####################################################################################################
# initializing the tssExp object                                                                   
initializeExp("Example RAMPAGE experiment (Rice)", "riceRAMPAGE", "bam_filtered", isPairedEnd=TRUE)#
####################################################################################################
# setting the sample IDs and names for the datasets within the RAMPAGE library
setSampleID(riceRAMPAGE, c("DL1_1","DL1_2", "DL2_1", "DL2_2", "DL2_3"), c(1,1,2,2,2))               #            #####################################################################################################
# loading the BAM files from the TSS profiling experiment (in this case, RAMPAGE)                     
importBam(riceRAMPAGE)                                                                             #
####################################################################################################
# converting BAM data into TSS information and attaching it to your tssExp object                  
bamToTSS(riceRAMPAGE)                                                                              #
####################################################################################################
# constructing a tss abundance matrix for a given dataset
tssExpr(expName=riceRAMPAGE, tssNum=1, writeTable=FALSE)                                           #
####################################################################################################
# constructing a tss abundance matrix for a given dataset
tssExpr(expName=riceRAMPAGE, tssNum=2, writeTable=FALSE)                                           #
####################################################################################################
# constructing a tss abundance matrix for a given dataset
tssExpr(expName=riceRAMPAGE, tssNum=3, writeTable=FALSE)                                           #
####################################################################################################
# constructing a tss abundance matrix for a given dataset
tssExpr(expName=riceRAMPAGE, tssNum=4, writeTable=FALSE)                                           #
####################################################################################################
# constructing a tss abundance matrix for a given dataset
tssExpr(expName=riceRAMPAGE, tssNum=5, writeTable=FALSE)                                           #
####################################################################################################
# finding TSRs for a given dataset (in this case, the first one)
tsrFind(expName=riceRAMPAGE, tssNum=1, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)	 
####################################################################################################
# finding TSRs for a given dataset (in this case, the second one)
tsrFind(expName=riceRAMPAGE, tssNum=2, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)       
####################################################################################################
# finding TSRs for a given dataset (in this case, the third one)
tsrFind(expName=riceRAMPAGE, tssNum=3, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)       
####################################################################################################
# finding TSRs for a given dataset (in this case, the fourth one)
tsrFind(expName=riceRAMPAGE, tssNum=4, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)       ####################################################################################################
# finding TSRs for a given dataset (in this case, the fifth one)
tsrFind(expName=riceRAMPAGE, tssNum=5, nTSSs=3, clustDist=20, setToCluster="replicates", writeTable=FALSE)	 
####################################################################################################
# merging tssData objects according to the replicate info
mergeTSS(expName=riceRAMPAGE)                                                            	   #
####################################################################################################
# merging exprData objects according to the replicate info
mergeExpr(expName=riceRAMPAGE)                                                            	   #
####################################################################################################
# finding TSRs from the merged datasets
tsrFind(expName=riceRAMPAGE, tssNum=1, nTSSs=3, clustDist=20, setToCluster="merged", writeTable=FALSE)	         
####################################################################################################
save(riceRAMPAGE, file="demo-RAMPAGE-test.RData") #uncomment if you'd like to save a binary of the tssExp object to your working directory
