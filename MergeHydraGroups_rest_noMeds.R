#################
### LOAD DATA ###
#################
subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/n763_DepAnxTd_rest_noMeds_subjData.rds")
hydraData <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/subjectData/HydraOutput/DepAnxTd_ctVolJLF_HydraSubtypes.csv", header=TRUE)

##################
### MERGE DATA ###
##################
dataComb <-merge(subjData, hydraData, by="bblid", all=FALSE)

n <- nrow(dataComb)

########################
### RECODE VARIABLES ###
########################

#Cluster membership: -1 denotes controls, positive integers denote subgroups. We need to change to 0=controls (comparison group in lm analyses), 1=subtype1, 2=subtype2, etc.
dataComb$Hydra_k1[which(dataComb$Hydra_k1==-1)] <- 0
dataComb$Hydra_k2[which(dataComb$Hydra_k2==-1)] <- 0
dataComb$Hydra_k3[which(dataComb$Hydra_k3==-1)] <- 0
dataComb$Hydra_k4[which(dataComb$Hydra_k4==-1)] <- 0
dataComb$Hydra_k5[which(dataComb$Hydra_k5==-1)] <- 0
dataComb$Hydra_k6[which(dataComb$Hydra_k6==-1)] <- 0
dataComb$Hydra_k7[which(dataComb$Hydra_k7==-1)] <- 0
dataComb$Hydra_k8[which(dataComb$Hydra_k8==-1)] <- 0
dataComb$Hydra_k9[which(dataComb$Hydra_k9==-1)] <- 0
dataComb$Hydra_k10[which(dataComb$Hydra_k10==-1)] <- 0

#Make the hydra group variables into factors. TD is the comparison group.
dataComb$Hydra_k1 <- as.factor(dataComb$Hydra_k1)
dataComb$Hydra_k2 <- as.factor(dataComb$Hydra_k2)
dataComb$Hydra_k3 <- as.factor(dataComb$Hydra_k3)
dataComb$Hydra_k4 <- as.factor(dataComb$Hydra_k4)
dataComb$Hydra_k5 <- as.factor(dataComb$Hydra_k5)
dataComb$Hydra_k6 <- as.factor(dataComb$Hydra_k6)
dataComb$Hydra_k7 <- as.factor(dataComb$Hydra_k7)
dataComb$Hydra_k8 <- as.factor(dataComb$Hydra_k8)
dataComb$Hydra_k9 <- as.factor(dataComb$Hydra_k9)
dataComb$Hydra_k10 <- as.factor(dataComb$Hydra_k10)

#Reorder the group variable of interest to make S2 the comparison group.
dataComb$Hydra_k2_reordered <- factor(dataComb$Hydra_k2, levels=c("2","1","0"))

#################
### SAVE DATA ###
#################
saveRDS(dataComb, "/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n763_DepAnxTd_rest_hydra_noMeds_subjData.rds")
