########################################
#### GROUP ANALYSES WITH DTI SAMPLE ####
########################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n923_DepAnxTd_dti_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)
library(schoolmath)

#Total sample size
n <- nrow(subjData)


########################
### Rename variables ###
########################

#The tracts are named slightly differently in the datafreeze dti variables than the JHU atlas names. Rename to match.
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_atr_l'] <- 'dti_dtitk_jhutract_fa_atr_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_atr_r'] <- 'dti_dtitk_jhutract_fa_atr_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cst_l'] <- 'dti_dtitk_jhutract_fa_cst_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cst_r'] <- 'dti_dtitk_jhutract_fa_cst_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cgc_l'] <- 'dti_dtitk_jhutract_fa_cgc_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cgc_r'] <- 'dti_dtitk_jhutract_fa_cgc_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cgh_l'] <- 'dti_dtitk_jhutract_fa_cgh_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_cgh_r'] <- 'dti_dtitk_jhutract_fa_cgh_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_ifo_l'] <- 'dti_dtitk_jhutract_fa_ifo_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_ifo_r'] <- 'dti_dtitk_jhutract_fa_ifo_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_ilf_l'] <- 'dti_dtitk_jhutract_fa_ilf_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_ilf_r'] <- 'dti_dtitk_jhutract_fa_ilf_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_slf_l'] <- 'dti_dtitk_jhutract_fa_slf_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_slf_r'] <- 'dti_dtitk_jhutract_fa_slf_R'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_uf_l'] <- 'dti_dtitk_jhutract_fa_uf_L'
names(subjData)[names(subjData) == 'dti_dtitk_jhutract_fa_uf_r'] <- 'dti_dtitk_jhutract_fa_uf_R'

##################
### DTI TRACTS ###
##################

#Get variable names
dtiTrRegions <- names(subjData)[grep("dti_dtitk_jhutract",names(subjData))]

##Run models (all three groups included). 
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
dtiTrGam <- lapply(dtiTrRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + dti64Tsnr + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
dtiTrGamSumm <- lapply(dtiTrGam, summary)


##Pass the gam results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
dtiTrAnova <- lapply(dtiTrGam, anova)

#Pull the p-values
dtiTr_p <- sapply(dtiTrAnova, function(v) v$pTerms.pv[3])

#Convert to data frame
dtiTr_p <- as.data.frame(dtiTr_p)

#Add row names for ease of viewing
rownames(dtiTr_p) <- dtiTrRegions

#Print original p-values to three decimal places
dtiTr_p_round <- round(dtiTr_p,3)


##Pull the F-values for the table
dtiTr_F <- sapply(dtiTrAnova, function(v) v$pTerms.table[3,2])

#Convert to data frame
dtiTr_F <- as.data.frame(dtiTr_F)

#Print F-values to two decimal places
dtiTr_F_round <- round(dtiTr_F,2)

#Add tract names
dtiTrNames <- as.data.frame(dtiTrRegions)
dtiTr_F_names <- cbind(dtiTrNames,dtiTr_F_round)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
dtiTr_p_fdr <- p.adjust(dtiTr_p[,1],method="fdr")

#Convert to data frame
dtiTr_p_fdr <- as.data.frame(dtiTr_p_fdr)

#Print fdr-corrected p-values to three decimal places
dtiTr_p_fdr_round <- round(dtiTr_p_fdr,3)

#Add region names
dtiTr_p_names <- cbind(dtiTrNames,dtiTr_p_fdr_round)

#Merge F values
dtiTr_omnibus <- merge(dtiTr_F_names,dtiTr_p_names, by="dtiTrRegions")

#Keep the original order of the regions
dtiTr_omnibus <- dtiTr_omnibus[dtiTr_p_names$dtiTrRegions,]

#Trim table to only regions that passed FDR correction
dtiTr_omnibus_signif <- dtiTr_omnibus[dtiTr_p_fdr<0.05,]

#Save omnibus results as a .csv file
write.table(dtiTr_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_dtiTr.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote=FALSE)


##Run same model with Hydra_k2_reordered to make S1 the comparison group in order to see S1 vs S2 differences in the gam model.
#This model will provide estimates and p-values for S1 vs. S2.
dtiTrGam_reordered <- lapply(dtiTrRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + dti64Tsnr + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
dtiTrGamSumm_reordered <- lapply(dtiTrGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
dtiTr_S1vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[4,4])
dtiTr_S2vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[5,4])
dtiTr_S1vsS2 <- sapply(dtiTrGam_reordered, function(v) summary(v)$p.table[4,4])

#Combine the pairwise p values
dtiTr_pairs <- cbind(dtiTr_S1vsTd,dtiTr_S2vsTd,dtiTr_S1vsS2)

#Convert to data frame
dtiTr_pairs <- as.data.frame(dtiTr_pairs)

#Add row names for ease of viewing
rownames(dtiTr_pairs) <- dtiTrRegions

#Print original p-values to three decimal places
dtiTr_pairs_round <- round(dtiTr_pairs,3)


##Pull the B-values for the table
dtiTr_B_S1vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[4,1])
dtiTr_B_S2vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[5,1])
dtiTr_B_S1vsS2 <- sapply(dtiTrGam_reordered, function(v) summary(v)$p.table[4,1])

#Convert to data frame
dtiTr_B_S1vsTd <- as.data.frame(dtiTr_B_S1vsTd)
dtiTr_B_S2vsTd <- as.data.frame(dtiTr_B_S2vsTd)
dtiTr_B_S1vsS2 <- as.data.frame(dtiTr_B_S1vsS2)

#Print B-values to two decimal places
dtiTr_B_S1vsTd_round <- round(dtiTr_B_S1vsTd,2)
dtiTr_B_S2vsTd_round <- round(dtiTr_B_S2vsTd,2)
dtiTr_B_S1vsS2_round <- round(dtiTr_B_S1vsS2,2)

#Add tract names
dtiTr_B_names <- cbind(dtiTrNames,dtiTr_B_S1vsTd_round,dtiTr_B_S2vsTd_round,dtiTr_B_S1vsS2_round)


##Pull the SE values for the table
dtiTr_SE_S1vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[4,2])
dtiTr_SE_S2vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[5,2])
dtiTr_SE_S1vsS2 <- sapply(dtiTrGam_reordered, function(v) summary(v)$p.table[4,2])

#Convert to data frame
dtiTr_SE_S1vsTd <- as.data.frame(dtiTr_SE_S1vsTd)
dtiTr_SE_S2vsTd <- as.data.frame(dtiTr_SE_S2vsTd)
dtiTr_SE_S1vsS2 <- as.data.frame(dtiTr_SE_S1vsS2)

#Print SE values to two decimal places
dtiTr_SE_S1vsTd_round <- round(dtiTr_SE_S1vsTd,2)
dtiTr_SE_S2vsTd_round <- round(dtiTr_SE_S2vsTd,2)
dtiTr_SE_S1vsS2_round <- round(dtiTr_SE_S1vsS2,2)

#Add tract names
dtiTr_SE_names <- cbind(dtiTrNames,dtiTr_SE_S1vsTd_round,dtiTr_SE_S2vsTd_round,dtiTr_SE_S1vsS2_round)


##Pull the t-values for the table
dtiTr_t_S1vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[4,3])
dtiTr_t_S2vsTd <- sapply(dtiTrGam, function(v) summary(v)$p.table[5,3])
dtiTr_t_S1vsS2 <- sapply(dtiTrGam_reordered, function(v) summary(v)$p.table[4,3])

#Convert to data frame
dtiTr_t_S1vsTd <- as.data.frame(dtiTr_t_S1vsTd)
dtiTr_t_S2vsTd <- as.data.frame(dtiTr_t_S2vsTd)
dtiTr_t_S1vsS2 <- as.data.frame(dtiTr_t_S1vsS2)

#Print t-values to two decimal places
dtiTr_t_S1vsTd_round <- round(dtiTr_t_S1vsTd,2)
dtiTr_t_S2vsTd_round <- round(dtiTr_t_S2vsTd,2)
dtiTr_t_S1vsS2_round <- round(dtiTr_t_S1vsS2,2)

#Add tract names
dtiTr_t_names <- cbind(dtiTrNames,dtiTr_t_S1vsTd_round,dtiTr_t_S2vsTd_round,dtiTr_t_S1vsS2_round)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
dtiTr_fdrTable<-as.data.frame(matrix(nrow=18,ncol=3))
colnames(dtiTr_fdrTable)[1]<-"dtiTr_S1vsTd_pfdr"
colnames(dtiTr_fdrTable)[2]<-"dtiTr_S2vsTd_pfdr"
colnames(dtiTr_fdrTable)[3]<-"dtiTr_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(dtiTr_pairs)) {
    row <- dtiTr_pairs[i,]
    dtiTr_fdrTable[i,] <- p.adjust(dtiTr_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
dtiTr_fdrTable_round <- round(dtiTr_fdrTable,3)

#Add region names
dtiTr_p_pairwise_names <- cbind(dtiTrNames,dtiTr_fdrTable_round)

#Merge B, SE, t, and pfdr values
dtiTr_B_SE <- merge(dtiTr_B_names,dtiTr_SE_names, by="dtiTrRegions")
dtiTr_B_SE_t <- merge(dtiTr_B_SE,dtiTr_t_names, by="dtiTrRegions")
dtiTr_B_SE_t_pfdr <- merge(dtiTr_B_SE_t,dtiTr_p_pairwise_names, by="dtiTrRegions")

#Keep the original order of the regions
dtiTr_pairwise <- dtiTr_B_SE_t_pfdr[dtiTr_p_pairwise_names$dtiTrRegions,]

#Reorder the columns
dtiTr_pairwise <- dtiTr_pairwise[,c(1,2,5,8,11,3,6,9,12,4,7,10,13)]

#Trim table to only regions that passed FDR correction for the omnibus test
dtiTr_pairwise_signif <- dtiTr_pairwise[dtiTr_p_fdr<0.05,]

#Save the pairwise results as a .csv file
write.table(dtiTr_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_dtiTr.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote=FALSE)


##Create table of tract names and t-values for mapping in fslview
#Load JHU labels and corresponding index numbers
#These are taken from /home/arosen/template/jlf/hiLoLookup/jhuTract.csv but Adon confirmed that tracts 19 and 20 were removed from the datafreeze variables because they were small and unreliable.

JHU.labels <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/images/jhuTract.csv", header=TRUE)

#Trim table to only regions that passed FDR correction for the omnibus test
dtiTr_t_S1vsS2_omnibus <- dtiTr_t_names[dtiTr_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs S2 pairwise test
dtiTr_t_S1vsS2_signif <- dtiTr_t_S1vsS2_omnibus[dtiTr_pairwise_signif$dtiTr_S1vsS2_pfdr<0.05,]

#Direction of results
dtiTr_t_S1vsS2_neg <- is.negative(dtiTr_t_S1vsS2_signif$dtiTr_t_S1vsS2)

#Rename variables for merging
dtiTr_t_S1vsS2_renamed <- rename(dtiTr_t_S1vsS2_signif, TRACT_NAME = dtiTrRegions)

#Keep only the S1vsS2 t-values
dtiTr_t_S1vsS2_renamed$dtiTr_t_S1vsTd <- NULL
dtiTr_t_S1vsS2_renamed$dtiTr_t_S2vsTd <- NULL

#Remove the "dti_dtitk_jhutract_fa_" prefix to match with the JHU labels
dtiTr_t_S1vsS2_renamed[] <- lapply(dtiTr_t_S1vsS2_renamed, function(x) gsub("dti_dtitk_jhutract_fa_", "", x))

#Merge to add index numbers that correspond to the significant tracts
dtiTr_merged_S1vsS2 <-merge(dtiTr_t_S1vsS2_renamed,JHU.labels, by="TRACT_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
dtiTr_merged_S1vsS2$TRACT_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only	makes sense if all your results are in the same direction (all negative).
dtiTr_merged_S1vsS2$dtiTr_t_S1vsS2 <- as.numeric(dtiTr_merged_S1vsS2$dtiTr_t_S1vsS2)
dtiTr_merged_S1vsS2$dtiTr_t_S1vsS2 <- abs(dtiTr_merged_S1vsS2$dtiTr_t_S1vsS2)

#Save as a .csv
write.table(dtiTr_merged_S1vsS2, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_dtiTr_S1vsS2.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Create table of tract names and t-values for mapping in fslview
#Trim table to only regions that passed FDR correction for the omnibus test
dtiTr_t_S1vsTd_omnibus <- dtiTr_t_names[dtiTr_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs TD pairwise test
dtiTr_t_S1vsTd_signif <- dtiTr_t_S1vsTd_omnibus[dtiTr_pairwise_signif$dtiTr_S1vsTd_pfdr<0.05,]

#Direction of results
dtiTr_t_S1vsTd_neg <- is.negative(dtiTr_t_S1vsTd_signif$dtiTr_t_S1vsTd)

#Rename variables for merging
dtiTr_t_S1vsTd_renamed <- rename(dtiTr_t_S1vsTd_signif, TRACT_NAME = dtiTrRegions)

#Keep only the S1vsTd t-values
dtiTr_t_S1vsTd_renamed$dtiTr_t_S2vsTd <- NULL
dtiTr_t_S1vsTd_renamed$dtiTr_t_S1vsS2 <- NULL

#Remove the "dti_dtitk_jhutract_fa_" prefix to match with the JHU labels
dtiTr_t_S1vsTd_renamed[] <- lapply(dtiTr_t_S1vsTd_renamed, function(x) gsub("dti_dtitk_jhutract_fa_", "", x))

#Merge to add index numbers that correspond to the significant tracts
dtiTr_merged_S1vsTd <-merge(dtiTr_t_S1vsTd_renamed,JHU.labels, by="TRACT_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
dtiTr_merged_S1vsTd$TRACT_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
dtiTr_merged_S1vsTd$dtiTr_t_S1vsTd <- as.numeric(dtiTr_merged_S1vsTd$dtiTr_t_S1vsTd)
dtiTr_merged_S1vsTd$dtiTr_t_S1vsTd <- abs(dtiTr_merged_S1vsTd$dtiTr_t_S1vsTd)

#Save as a .csv
write.table(dtiTr_merged_S1vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_dtiTr_S1vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Create table of tract names and t-values for mapping in fslview
#Trim table to only regions that passed FDR correction for the omnibus test
dtiTr_t_S2vsTd_omnibus <- dtiTr_t_names[dtiTr_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S2 vs TD pairwise test
dtiTr_t_S2vsTd_signif <- dtiTr_t_S2vsTd_omnibus[dtiTr_pairwise_signif$dtiTr_S2vsTd_pfdr<0.05,]

#Direction of results
dtiTr_t_S2vsTd_neg <- is.negative(dtiTr_t_S2vsTd_signif$dtiTr_t_S2vsTd)

#Rename variables for merging
dtiTr_t_S2vsTd_renamed <- rename(dtiTr_t_S2vsTd_signif, TRACT_NAME = dtiTrRegions)

#Keep only the S2vsTd t-values
dtiTr_t_S2vsTd_renamed$dtiTr_t_S1vsTd <- NULL
dtiTr_t_S2vsTd_renamed$dtiTr_t_S1vsS2 <- NULL

#Remove the "dti_dtitk_jhutract_fa_" prefix to match with the JHU labels
dtiTr_t_S2vsTd_renamed[] <- lapply(dtiTr_t_S2vsTd_renamed, function(x) gsub("dti_dtitk_jhutract_fa_", "", x))

#Merge to add index numbers that correspond to the significant tracts
dtiTr_merged_S2vsTd <-merge(dtiTr_t_S2vsTd_renamed,JHU.labels, by="TRACT_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
dtiTr_merged_S2vsTd$TRACT_NAME <- NULL

#Save as a .csv
write.table(dtiTr_merged_S2vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_dtiTr_S2vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")
