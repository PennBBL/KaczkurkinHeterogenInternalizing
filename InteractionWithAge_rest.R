#######################################################
#### INTERACTIONS WITH AGE IN RESTING STATE SAMPLE ####
#######################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n840_DepAnxTd_rest_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)

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
  gam(substitute(i ~ s(age) + sex + restRelMeanRMSMotion + Hydra_k2 + Hydra_k2*age, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
restAlffGamSumm <- lapply(restAlffGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
restAlffAnova <- lapply(restAlffGam, anova)

#Pull the p-values
restAlff_p <- sapply(restAlffAnova, function(v) v$pTerms.pv[5])

#Convert to data frame
restAlff_p <- as.data.frame(restAlff_p)

#Add row names for ease of viewing
rownames(restAlff_p) <- restAlffRegions

#Print original p-values to three decimal places
restAlff_p_round <- round(restAlff_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
restAlff_p_fdr <- p.adjust(restAlff_p[,1],method="fdr")

#Convert to data frame
restAlff_p_fdr <- as.data.frame(restAlff_p_fdr)

#Print fdr-corrected p-values to three decimal places
restAlff_p_fdr_round <- round(restAlff_p_fdr,3)

#Add region names
restAlffNames <- as.data.frame(restAlffRegions)
restAlff_p_names <- cbind(restAlffNames,restAlff_p_fdr_round)

#Trim table to only regions that passed FDR correction for the omnibus test
restAlff_signif <- restAlff_p_names[restAlff_p_fdr<0.05,]
