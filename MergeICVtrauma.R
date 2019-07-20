#################
### LOAD DATA ###
#################
subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_subjData.rds")
icvData <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/subjectData/GO-BBL_muse_dramms+ants_C1.2_Features_DerivedVolumes.csv", header=TRUE)
traumaData <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/subjectData/traumaScores.csv", header=TRUE)

##################
### MERGE DATA ###
##################
#Pull out ICV
icv <- icvData[c(grep("ID|^X702$",names(icvData)))]

#ID contains bblid_scanid. Separate the two.
icv$bblid <- sub("_.*", "", icv$ID)
icv$scanid <- sub(".*?_", "", icv$ID)

#icvData also has leading 0s in front of bblids and scanids (e.g., "080010"), while subjData does not ("80010"). Remove leading 0s.
icv$bblid <- as.numeric(icv$bblid)
icv$scanid <- as.numeric(icv$scanid)

#Rename variables
names(icv)[names(icv) == 'X702'] <- 'ICV'

#Remove ID
icv$ID <- NULL

#Merge but keep only 1141 sample
dataComb <-merge(subjData, icv, by=c("bblid","scanid"), all=FALSE)

#Merge trauma scores but keep only the 1141 sample
dataComb2 <-merge(dataComb, traumaData, by="bblid", all=FALSE)

#Put bblids in ascending order
data.ordered <- dataComb2[order(dataComb2$bblid),]

n <- nrow(data.ordered)

#################
### SAVE DATA ###
#################
saveRDS(data.ordered, "/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_ICV_trauma_subjData.rds")
