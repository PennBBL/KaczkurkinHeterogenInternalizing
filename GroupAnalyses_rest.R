##################################################
#### GROUP ANALYSES WITH RESTING STATE SAMPLE ####
##################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n840_DepAnxTd_rest_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)
library(schoolmath)

#Total sample size
n <- nrow(subjData)


######################
### Rest ALFF ROIs ###
######################

#Get variable names
restAlff_all <- subjData[c(grep("rest_jlf_alff",names(subjData)))]
restAlff_short <- restAlff_all[,-grep("Cerebell",names(restAlff_all))]
restAlffRegions <- names(restAlff_short)

##Run models (all three groups included). 
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
restAlffGam <- lapply(restAlffRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + restRelMeanRMSMotion + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
restAlffGamSumm <- lapply(restAlffGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
restAlffAnova <- lapply(restAlffGam, anova)

#Pull the p-values
restAlff_p <- sapply(restAlffAnova, function(v) v$pTerms.pv[3])

#Convert to data frame
restAlff_p <- as.data.frame(restAlff_p)

#Add row names for ease of viewing
rownames(restAlff_p) <- restAlffRegions

#Print original p-values to three decimal places
restAlff_p_round <- round(restAlff_p,3)


##Pull the F-values for the table
restAlff_F <- sapply(restAlffAnova, function(v) v$pTerms.table[3,2])

#Convert to data frame
restAlff_F <- as.data.frame(restAlff_F)

#Print F-values to two decimal places
restAlff_F_round <- round(restAlff_F,2)

#Add region names
restAlffNames <- as.data.frame(restAlffRegions)
restAlff_F_names <- cbind(restAlffNames,restAlff_F_round)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
restAlff_p_fdr <- p.adjust(restAlff_p[,1],method="fdr")

#Convert to data frame
restAlff_p_fdr <- as.data.frame(restAlff_p_fdr)

#Print fdr-corrected p-values to three decimal places
restAlff_p_fdr_round <- round(restAlff_p_fdr,3)

#Add region names
restAlff_p_names <- cbind(restAlffNames,restAlff_p_fdr_round)

#Merge F values
restAlff_omnibus <- merge(restAlff_F_names,restAlff_p_names, by="restAlffRegions")

#Keep the original order of the regions
restAlff_omnibus <- restAlff_omnibus[restAlff_p_names$restAlffRegions,]

#Trim table to only regions that passed FDR correction
restAlff_omnibus_signif <- restAlff_omnibus[restAlff_p_fdr<0.05,]

#Save omnibus results as a .csv file
write.table(restAlff_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_restAlff.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
restAlffGam_reordered <- lapply(restAlffRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + restRelMeanRMSMotion + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
restAlffGamSumm_reordered <- lapply(restAlffGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
restAlff_S1vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[4,4])
restAlff_S2vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[5,4])
restAlff_S1vsS2 <- sapply(restAlffGam_reordered, function(v) summary(v)$p.table[4,4])

#Combine the pairwise p values
restAlff_pairs <- cbind(restAlff_S1vsTd,restAlff_S2vsTd,restAlff_S1vsS2)

#Convert to data frame
restAlff_pairs <- as.data.frame(restAlff_pairs)

#Add row names for ease of viewing
rownames(restAlff_pairs) <- restAlffRegions

#Print original p-values to three decimal places
restAlff_pairs_round <- round(restAlff_pairs,3)


##Pull the B-values for the table
restAlff_B_S1vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[4,1])
restAlff_B_S2vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[5,1])
restAlff_B_S1vsS2 <- sapply(restAlffGam_reordered, function(v) summary(v)$p.table[4,1])

#Convert to data frame
restAlff_B_S1vsTd <- as.data.frame(restAlff_B_S1vsTd)
restAlff_B_S2vsTd <- as.data.frame(restAlff_B_S2vsTd)
restAlff_B_S1vsS2 <- as.data.frame(restAlff_B_S1vsS2)

#Print B-values to two decimal places
restAlff_B_S1vsTd_round <- round(restAlff_B_S1vsTd,2)
restAlff_B_S2vsTd_round <- round(restAlff_B_S2vsTd,2)
restAlff_B_S1vsS2_round <- round(restAlff_B_S1vsS2,2)

#Add tract names
restAlff_B_names <- cbind(restAlffNames,restAlff_B_S1vsTd_round,restAlff_B_S2vsTd_round,restAlff_B_S1vsS2_round)


##Pull the SE values for the table
restAlff_SE_S1vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[4,2])
restAlff_SE_S2vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[5,2])
restAlff_SE_S1vsS2 <- sapply(restAlffGam_reordered, function(v) summary(v)$p.table[4,2])

#Convert to data frame
restAlff_SE_S1vsTd <- as.data.frame(restAlff_SE_S1vsTd)
restAlff_SE_S2vsTd <- as.data.frame(restAlff_SE_S2vsTd)
restAlff_SE_S1vsS2 <- as.data.frame(restAlff_SE_S1vsS2)

#Print SE values to two decimal places
restAlff_SE_S1vsTd_round <- round(restAlff_SE_S1vsTd,2)
restAlff_SE_S2vsTd_round <- round(restAlff_SE_S2vsTd,2)
restAlff_SE_S1vsS2_round <- round(restAlff_SE_S1vsS2,2)

#Add tract names
restAlff_SE_names <- cbind(restAlffNames,restAlff_SE_S1vsTd_round,restAlff_SE_S2vsTd_round,restAlff_SE_S1vsS2_round)


##Pull the t-values for the table
restAlff_t_S1vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[4,3])
restAlff_t_S2vsTd <- sapply(restAlffGam, function(v) summary(v)$p.table[5,3])
restAlff_t_S1vsS2 <- sapply(restAlffGam_reordered, function(v) summary(v)$p.table[4,3])

#Convert to data frame
restAlff_t_S1vsTd <- as.data.frame(restAlff_t_S1vsTd)
restAlff_t_S2vsTd <- as.data.frame(restAlff_t_S2vsTd)
restAlff_t_S1vsS2 <- as.data.frame(restAlff_t_S1vsS2)

#Print t-values to two decimal places
restAlff_t_S1vsTd_round <- round(restAlff_t_S1vsTd,2)
restAlff_t_S2vsTd_round <- round(restAlff_t_S2vsTd,2)
restAlff_t_S1vsS2_round <- round(restAlff_t_S1vsS2,2)

#Add tract names
restAlff_t_names <- cbind(restAlffNames,restAlff_t_S1vsTd_round,restAlff_t_S2vsTd_round,restAlff_t_S1vsS2_round)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
restAlff_fdrTable<-as.data.frame(matrix(nrow=112,ncol=3))
colnames(restAlff_fdrTable)[1]<-"restAlff_S1vsTd_pfdr"
colnames(restAlff_fdrTable)[2]<-"restAlff_S2vsTd_pfdr"
colnames(restAlff_fdrTable)[3]<-"restAlff_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(restAlff_pairs)) {
    row <- restAlff_pairs[i,]
    restAlff_fdrTable[i,] <- p.adjust(restAlff_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
restAlff_fdrTable_round <- round(restAlff_fdrTable,3)

#Add region names
restAlff_p_pairwise_names <- cbind(restAlffNames,restAlff_fdrTable_round)

#Merge B, SE, t, and pfdr values
restAlff_B_SE <- merge(restAlff_B_names,restAlff_SE_names, by="restAlffRegions")
restAlff_B_SE_t <- merge(restAlff_B_SE,restAlff_t_names, by="restAlffRegions")
restAlff_B_SE_t_pfdr <- merge(restAlff_B_SE_t,restAlff_p_pairwise_names, by="restAlffRegions")

#Keep the original order of the regions
restAlff_pairwise <- restAlff_B_SE_t_pfdr[restAlff_p_pairwise_names$restAlffRegions,]

#Reorder the columns
restAlff_pairwise <- restAlff_pairwise[,c(1,2,5,8,11,3,6,9,12,4,7,10,13)]

#Trim table to only regions that passed FDR correction for the omnibus test
restAlff_pairwise_signif <- restAlff_pairwise[restAlff_p_fdr<0.05,]

#Save the pairwise results as a .csv file
write.table(restAlff_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_restAlff.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Create table of region names and t-values for mapping in fslview
#Load JLF labels and corresponding index numbers
#These are taken from https://github.com/PennBBL/jlfVisualizer/blob/master/data/jlf_lookupWithWM.csv

JLF.labels <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/images/jlf_lookupWithWM.csv", header=TRUE)

#Trim table to only regions that passed FDR correction for the omnibus test
restAlff_t_S1vsS2_omnibus <- restAlff_t_names[restAlff_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs S2 pairwise test
restAlff_t_S1vsS2_signif <- restAlff_t_S1vsS2_omnibus[restAlff_pairwise_signif$restAlff_S1vsS2_pfdr<0.05,]

#Direction of results
restAlff_t_S1vsS2_neg <- is.negative(restAlff_t_S1vsS2_signif$restAlff_t_S1vsS2)

#Rename variables for merging
restAlff_t_S1vsS2_renamed <- rename(restAlff_t_S1vsS2_signif, ROI_NAME = restAlffRegions)

#Keep only the S1vsS2 t-values
restAlff_t_S1vsS2_renamed$restAlff_t_S1vsTd <- NULL
restAlff_t_S1vsS2_renamed$restAlff_t_S2vsTd <- NULL

#Remove the "rest_jlf_alff_" prefix to match with the JLF labels
restAlff_t_S1vsS2_renamed[] <- lapply(restAlff_t_S1vsS2_renamed, function(x) gsub("rest_jlf_alff_", "", x))

#Merge to add index numbers that correspond to the significant regions
restAlff_merged_S1vsS2 <-merge(restAlff_t_S1vsS2_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
restAlff_merged_S1vsS2$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only	makes sense if all your results are in the same direction (all negative).
restAlff_merged_S1vsS2$restAlff_t_S1vsS2 <- as.numeric(restAlff_merged_S1vsS2$restAlff_t_S1vsS2)
restAlff_merged_S1vsS2$restAlff_t_S1vsS2 <- abs(restAlff_merged_S1vsS2$restAlff_t_S1vsS2)

#Save as a .csv
write.table(restAlff_merged_S1vsS2, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_restAlff_S1vsS2.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Create table of region names and t-values for mapping in fslview
#Trim table to only regions that passed FDR correction for the omnibus test
restAlff_t_S1vsTd_omnibus <- restAlff_t_names[restAlff_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs TD pairwise test
restAlff_t_S1vsTd_signif <- restAlff_t_S1vsTd_omnibus[restAlff_pairwise_signif$restAlff_S1vsTd_pfdr<0.05,]

#Direction of results
restAlff_t_S1vsTd_neg <- is.negative(restAlff_t_S1vsTd_signif$restAlff_t_S1vsTd)

#Rename variables for merging
restAlff_t_S1vsTd_renamed <- rename(restAlff_t_S1vsTd_signif, ROI_NAME = restAlffRegions)

#Keep only the S1vsTd t-values
restAlff_t_S1vsTd_renamed$restAlff_t_S2vsTd <- NULL
restAlff_t_S1vsTd_renamed$restAlff_t_S1vsS2 <- NULL

#Remove the "rest_jlf_alff_" prefix to match with the JLF labels
restAlff_t_S1vsTd_renamed[] <- lapply(restAlff_t_S1vsTd_renamed, function(x) gsub("rest_jlf_alff_", "", x))

#Merge to add index numbers that correspond to the significant regions
restAlff_merged_S1vsTd <-merge(restAlff_t_S1vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
restAlff_merged_S1vsTd$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
restAlff_merged_S1vsTd$restAlff_t_S1vsTd <- as.numeric(restAlff_merged_S1vsTd$restAlff_t_S1vsTd)
restAlff_merged_S1vsTd$restAlff_t_S1vsTd <- abs(restAlff_merged_S1vsTd$restAlff_t_S1vsTd)

#Save as a .csv
write.table(restAlff_merged_S1vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_restAlff_S1vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Create table of region names and t-values for mapping in fslview
#Trim table to only regions that passed FDR correction for the omnibus test
restAlff_t_S2vsTd_omnibus <- restAlff_t_names[restAlff_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S2 vs TD pairwise test
restAlff_t_S2vsTd_signif <- restAlff_t_S2vsTd_omnibus[restAlff_pairwise_signif$restAlff_S2vsTd_pfdr<0.05,]

#Direction of results
restAlff_t_S2vsTd_neg <- is.negative(restAlff_t_S2vsTd_signif$restAlff_t_S2vsTd)

#Rename variables for merging
restAlff_t_S2vsTd_renamed <- rename(restAlff_t_S2vsTd_signif, ROI_NAME = restAlffRegions)

#Keep only the S2vsTd t-values
restAlff_t_S2vsTd_renamed$restAlff_t_S1vsTd <- NULL
restAlff_t_S2vsTd_renamed$restAlff_t_S1vsS2 <- NULL

#Remove the "rest_jlf_alff_" prefix to match with the JLF labels
restAlff_t_S2vsTd_renamed[] <- lapply(restAlff_t_S2vsTd_renamed, function(x) gsub("rest_jlf_alff_", "", x))

#Merge to add index numbers that correspond to the significant regions
restAlff_merged_S2vsTd <-merge(restAlff_t_S2vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
restAlff_merged_S2vsTd$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
restAlff_merged_S2vsTd$restAlff_t_S2vsTd <- as.numeric(restAlff_merged_S2vsTd$restAlff_t_S2vsTd)
restAlff_merged_S2vsTd$restAlff_t_S2vsTd <- abs(restAlff_merged_S2vsTd$restAlff_t_S2vsTd)

#Save as a .csv
write.table(restAlff_merged_S2vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_restAlff_S2vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")
