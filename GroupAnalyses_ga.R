################################
#### GROUP ANALYSES WITH GA ####
################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n232_DepAnxTd_t1_hydra_GA_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)
library(psych)
library(schoolmath)

#Total sample size
n <- nrow(subjData)

#######################
### Gestational Age ###
#######################

##Run model
gaLm <- lm(ga ~ Hydra_k2, data=subjData)

#Look at model summary
gaLmSumm <- summary(gaLm)


##Pass the results to the anova() function to fit ANOVAs
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
gaAnova <- anova(gaLm)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
gaLm_reordered <- lm(ga ~ Hydra_k2_reordered, data=subjData)

#Look at model summary
gaLmSumm_reordered <- summary(gaLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
ga_S1vsTd <- summary(gaLm)$coefficients[2,4]
ga_S2vsTd <- summary(gaLm)$coefficients[3,4]
ga_S1vsS2 <- summary(gaLm_reordered)$coefficients[2,4]

#Combine the pairwise p values
ga_pairs <- cbind(ga_S1vsTd,ga_S2vsTd,ga_S1vsS2)

#Convert to data frame
ga_pairs <- as.data.frame(ga_pairs)

#Print original p-values to three decimal places
ga_pairs_round <- round(ga_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
ga_fdr <- p.adjust(ga_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
ga_fdr_round <- round(ga_fdr,3)


##Follow-up on significant results
#Means and SDs
gaDescriptives <- tapply(subjData$ga, subjData$Hydra_k2, describe)
group_n <- table(subjData$Hydra_k2)
