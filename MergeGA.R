#################
### LOAD DATA ###
#################
subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_ICV_trauma_subjData.rds")
gaData <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/subjectData/gaData_final.csv", header=TRUE)

##################
### MERGE DATA ###
##################
#Merge
dataComb <-merge(gaData, subjData, by="bblid", all=FALSE)

#Remove missing ga data
dataComb <- dataComb[!is.na(dataComb$ga),]

#Put bblids in ascending order
data.ordered <- dataComb[order(dataComb$bblid),]

n <- nrow(data.ordered)

#################
### SAVE DATA ###
#################
saveRDS(data.ordered, "/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n232_DepAnxTd_t1_hydra_GA_subjData.rds")
