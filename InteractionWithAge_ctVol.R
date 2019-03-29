#######################################################
#### INTERACTIONS WITH AGE IN VOLUME/CT JLF SAMPLE ####
#######################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)

#Total sample size
n <- nrow(subjData)


##############################
### JLF cortical thickness ###
##############################

#Get variable names
ctRegions <- names(subjData)[grep("mprage_jlf_ct",names(subjData))]

##Run models (all three groups included). 
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
ctGam <- lapply(ctRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2 + Hydra_k2*age, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
ctGamSumm <- lapply(ctGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
ctAnova <- lapply(ctGam, anova)

#Pull the p-values
ct_p <- sapply(ctAnova, function(v) v$pTerms.pv[5])

#Convert to data frame
ct_p <- as.data.frame(ct_p)

#Add row names for ease of viewing
rownames(ct_p) <- ctRegions

#Print original p-values to three decimal places
ct_p_round <- round(ct_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
ct_p_fdr <- p.adjust(ct_p[,1],method="fdr")

#Convert to data frame
ct_p_fdr <- as.data.frame(ct_p_fdr)

#Print fdr-corrected p-values to three decimal places
ct_p_fdr_round <- round(ct_p_fdr,3)

#Add region names
ctNames <- as.data.frame(ctRegions)
ct_p_names <- cbind(ctNames,ct_p_fdr_round)

#Trim table to only regions that passed FDR correction for the omnibus test
ct_signif <- ct_p_names[ct_p_fdr<0.05,]


##################
### JLF volume ###
##################

#Get variable names
volRegions_all <- subjData[c(grep("mprage_jlf_vol",names(subjData)))]
volRegions_short <- volRegions_all[,-grep("Vent|Brain_Stem|Cerebell|Cerebral_White_Matter|CSF|Lobe_WM",names(volRegions_all))]
volRegions <- names(volRegions_short)

##Run models (all three groups included).
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
volGam <- lapply(volRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2 + Hydra_k2*age, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
volGamSumm <- lapply(volGam, summary)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
volAnova <- lapply(volGam, anova)

#Pull the p-values
vol_p <- sapply(volAnova, function(v) v$pTerms.pv[5])

#Convert to data frame
vol_p <- as.data.frame(vol_p)

#Add row names for ease of viewing
rownames(vol_p) <- volRegions

#Print original p-values to three decimal places
vol_p_round <- round(vol_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values
vol_p_fdr <- p.adjust(vol_p[,1],method="fdr")

#Convert to data frame
vol_p_fdr <- as.data.frame(vol_p_fdr)

#Print fdr-corrected p-values to three decimal places
vol_p_fdr_round <- round(vol_p_fdr,3)

#Add region names
volNames <- as.data.frame(volRegions)
vol_p_names <- cbind(volNames,vol_p_fdr_round)

#Trim table to only regions that passed FDR correction for the omnibus test
vol_signif <- vol_p_names[vol_p_fdr<0.05,]
