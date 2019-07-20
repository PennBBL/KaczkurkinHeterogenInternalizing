################################
#### GROUP ANALYSES WITH DX ####
################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)
library(schoolmath)

#Total sample size
n <- nrow(subjData)


#############################################
### Psychopathology Diagnostic Categories ###
#############################################

#Get variable names
dx_all <- subjData[c(grep("goassessSmry",names(subjData)))]
dx_short <- dx_all[,-grep("Mood|Eat|Anx|Del|Hal|Beh|Prime|Overall",names(dx_all))]
dx <- names(dx_short)

##Run models (all three groups included). 
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
dxLm <- lapply(dx, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
dxLmSumm <- lapply(dxLm, summary)


##Pass the results to the anova() function to fit ANOVAs
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
dxAnova <- lapply(dxLm, anova)

#Pull the p-values
dx_p <- sapply(dxAnova, function(v) v$"Pr(>F)"[3])

#Convert to data frame
dx_p <- as.data.frame(dx_p)

#Add row names for ease of viewing
rownames(dx_p) <- dx

#Print original p-values to three decimal places
dx_p_round <- round(dx_p,3)


##Pull the F-values for the table
dx_F <- sapply(dxAnova, function(v) v$"F value"[3])

#Convert to data frame
dx_F <- as.data.frame(dx_F)

#Print F-values to two decimal places
dx_F_round <- round(dx_F,2)

#Add region names
dxNames <- as.data.frame(dx)
dx_F_names <- cbind(dxNames,dx_F_round)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
dx_p_fdr <- p.adjust(dx_p[,1],method="fdr")

#Convert to data frame
dx_p_fdr <- as.data.frame(dx_p_fdr)

#Print fdr-corrected p-values to three decimal places
dx_p_fdr_round <- round(dx_p_fdr,3)

#Add region names
dx_p_names <- cbind(dxNames,dx_p_fdr_round)

#Merge F values
dx_omnibus <- merge(dx_F_names,dx_p_names, by="dx")

#Keep the original order of the variables
dx_omnibus <- dx_omnibus[dx_p_names$dx,]

#Trim table to only variables that passed FDR correction
dx_omnibus_signif <- dx_omnibus[dx_p_fdr<0.05,]

#Save omnibus results as a .csv file
write.table(dx_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_dx.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
dxLm_reordered <- lapply(dx, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2_reordered, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
dxLmSumm_reordered <- lapply(dxLm_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
dx_S1vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[4,4])
dx_S2vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[5,4])
dx_S1vsS2 <- sapply(dxLm_reordered, function(v) summary(v)$coefficients[4,4])

#Combine the pairwise p values
dx_pairs <- cbind(dx_S1vsTd,dx_S2vsTd,dx_S1vsS2)

#Convert to data frame
dx_pairs <- as.data.frame(dx_pairs)

#Add row names for ease of viewing
rownames(dx_pairs) <- dx

#Print original p-values to three decimal places
dx_pairs_round <- round(dx_pairs,3)


##Pull the B-values for the table
dx_B_S1vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[4,1])
dx_B_S2vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[5,1])
dx_B_S1vsS2 <- sapply(dxLm_reordered, function(v) summary(v)$coefficients[4,1])

#Convert to data frame
dx_B_S1vsTd <- as.data.frame(dx_B_S1vsTd)
dx_B_S2vsTd <- as.data.frame(dx_B_S2vsTd)
dx_B_S1vsS2 <- as.data.frame(dx_B_S1vsS2)

#Print B-values to two decimal places
dx_B_S1vsTd_round <- round(dx_B_S1vsTd,2)
dx_B_S2vsTd_round <- round(dx_B_S2vsTd,2)
dx_B_S1vsS2_round <- round(dx_B_S1vsS2,2)

#Add names
dx_B_names <- cbind(dxNames,dx_B_S1vsTd_round,dx_B_S2vsTd_round,dx_B_S1vsS2_round)


##Pull the SE values for the table
dx_SE_S1vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[4,2])
dx_SE_S2vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[5,2])
dx_SE_S1vsS2 <- sapply(dxLm_reordered, function(v) summary(v)$coefficients[4,2])

#Convert to data frame
dx_SE_S1vsTd <- as.data.frame(dx_SE_S1vsTd)
dx_SE_S2vsTd <- as.data.frame(dx_SE_S2vsTd)
dx_SE_S1vsS2 <- as.data.frame(dx_SE_S1vsS2)

#Print SE values to two decimal places
dx_SE_S1vsTd_round <- round(dx_SE_S1vsTd,2)
dx_SE_S2vsTd_round <- round(dx_SE_S2vsTd,2)
dx_SE_S1vsS2_round <- round(dx_SE_S1vsS2,2)

#Add names
dx_SE_names <- cbind(dxNames,dx_SE_S1vsTd_round,dx_SE_S2vsTd_round,dx_SE_S1vsS2_round)


##Pull the t-values for the table
dx_t_S1vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[4,3])
dx_t_S2vsTd <- sapply(dxLm, function(v) summary(v)$coefficients[5,3])
dx_t_S1vsS2 <- sapply(dxLm_reordered, function(v) summary(v)$coefficients[4,3])

#Convert to data frame
dx_t_S1vsTd <- as.data.frame(dx_t_S1vsTd)
dx_t_S2vsTd <- as.data.frame(dx_t_S2vsTd)
dx_t_S1vsS2 <- as.data.frame(dx_t_S1vsS2)

#Print t-values to two decimal places
dx_t_S1vsTd_round <- round(dx_t_S1vsTd,2)
dx_t_S2vsTd_round <- round(dx_t_S2vsTd,2)
dx_t_S1vsS2_round <- round(dx_t_S1vsS2,2)

#Add names
dx_t_names <- cbind(dxNames,dx_t_S1vsTd_round,dx_t_S2vsTd_round,dx_t_S1vsS2_round)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
dx_fdrTable<-as.data.frame(matrix(nrow=16,ncol=3))
colnames(dx_fdrTable)[1]<-"dx_S1vsTd_pfdr"
colnames(dx_fdrTable)[2]<-"dx_S2vsTd_pfdr"
colnames(dx_fdrTable)[3]<-"dx_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(dx_pairs)) {
    row <- dx_pairs[i,]
    dx_fdrTable[i,] <- p.adjust(dx_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
dx_fdrTable_round <- round(dx_fdrTable,3)

#Add region names
dx_p_pairwise_names <- cbind(dxNames,dx_fdrTable_round)

#Merge B, SE, t, and pfdr values
dx_B_SE <- merge(dx_B_names,dx_SE_names, by="dx")
dx_B_SE_t <- merge(dx_B_SE,dx_t_names, by="dx")
dx_B_SE_t_pfdr <- merge(dx_B_SE_t,dx_p_pairwise_names, by="dx")

#Keep the original order of the variables
dx_pairwise <- dx_B_SE_t_pfdr[dx_p_pairwise_names$dx,]

#Reorder the columns
dx_pairwise <- dx_pairwise[,c(1,2,5,8,11,3,6,9,12,4,7,10,13)]

#Trim table to only variables that passed FDR correction for the omnibus test
dx_pairwise_signif <- dx_pairwise[dx_p_fdr<0.05,]

#Save the pairwise results as a .csv file
write.table(dx_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_dx.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)

