################################################
#### INTERACTIONS WITH GENDER IN DTI SAMPLE ####
################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n923_DepAnxTd_dti_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(dplyr)

#Total sample size
n <- nrow(subjData)


##################
### DTI TRACTS ###
##################

#Get variable names
dtiTrRegions <- names(subjData)[grep("dti_dtitk_jhutract",names(subjData))]

##Run models (all three groups included). 
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
dtiTrGam <- lapply(dtiTrRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + dti64Tsnr + Hydra_k2 + Hydra_k2*sex, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
dtiTrGamSumm <- lapply(dtiTrGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
dtiTrAnova <- lapply(dtiTrGam, anova)

#Pull the p-values
dtiTr_p <- sapply(dtiTrAnova, function(v) v$pTerms.pv[4])

#Convert to data frame
dtiTr_p <- as.data.frame(dtiTr_p)

#Add row names for ease of viewing
rownames(dtiTr_p) <- dtiTrRegions

#Print original p-values to three decimal places
dtiTr_p_round <- round(dtiTr_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values 
dtiTr_p_fdr <- p.adjust(dtiTr_p[,1],method="fdr")

#Convert to data frame
dtiTr_p_fdr <- as.data.frame(dtiTr_p_fdr)

#Print fdr-corrected p-values to three decimal places
dtiTr_p_fdr_round <- round(dtiTr_p_fdr,3)

#Add region names
dtiTrNames <- as.data.frame(dtiTrRegions)
dtiTr_p_names <- cbind(dtiTrNames,dtiTr_p_fdr_round)

#Trim table to only regions that passed FDR correction for the omnibus test
dtiTr_signif <- dtiTr_p_names[dtiTr_p_fdr<0.05,]
