##############################################
#### Prepare Subject-Level Data For HYDRA ####
##############################################

##Read in subject level data
DepAnxTd_t1 <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/n1141_DepAnxTd_t1_subjData.rds")

####################
#### Covariates ####
####################

#T1
DepAnxTd_t1_AgeSex <- DepAnxTd_t1[c(grep("bblid|^age$|sex",names(DepAnxTd_t1)))]

#######################################
#### Brain Features: CT and Volume ####
#######################################

#CT JLF
DepAnxTd_ctJLF <- DepAnxTd_t1[c(grep("bblid|mprage_jlf_ct",names(DepAnxTd_t1)))]

#Volume JLF (we removed ventricles, brain stem, cerebellum, white matter, CSF, and lobes)
DepAnxTd_volJLF_all <- DepAnxTd_t1[c(grep("bblid|mprage_jlf_vol|^DepAnxTd2$",names(DepAnxTd_t1)))]
DepAnxTd_volJLF <- DepAnxTd_volJLF_all[,-grep("Vent|Brain_Stem|Cerebell|Cerebral_White_Matter|CSF|Lobe_WM",names(DepAnxTd_volJLF_all))]

#CT and volume JLF
DepAnxTd_ctVolJLF <- merge(DepAnxTd_ctJLF,DepAnxTd_volJLF, all=TRUE)

#Put bblids in ascending order
DepAnxTd_ctVolJLF <- DepAnxTd_ctVolJLF[order(DepAnxTd_ctVolJLF$bblid),]

#####################
#### Check Files ####
#####################

##The covariates and features files MUST have subjects in the same order. Check that this true (correlation between bblids should = 1)
DepAnxTd_corr_ctVolJLF <- cor(DepAnxTd_t1_AgeSex$bblid,DepAnxTd_ctVolJLF$bblid)

##Check that the variables are in the correct order with the correct n
#n should be 1141
nrow(DepAnxTd_t1_AgeSex)
nrow(DepAnxTd_ctVolJLF)

#bblid needs to be the first variable; group needs to be the last variable
names(DepAnxTd_t1_AgeSex)
names(DepAnxTd_ctVolJLF)

####################
#### Save Files ####
####################

#Save covariates files to use with HYDRA
write.csv(DepAnxTd_t1_AgeSex, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/SubsetsForHYDRA/DepAnxTd_t1_AgeSex.csv", row.names=F, quote=F)

#Save brain features to use with HYDRA
write.csv(DepAnxTd_ctVolJLF, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/SubsetsForHYDRA/DepAnxTd_ctVolJLF.csv", row.names=F, quote=F)
