############################################################
#### GROUP ANALYSES WITH VOLUME/CT JLF SAMPLE WITH RACE ####
############################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_ICV_subjData.rds")

#Load libraries ("plyr" needs to be loaded before "dplyr")
library(plyr)
library(mgcv)
library(dplyr)
library(varhandle)
library(psych)
library(schoolmath)

#Total sample size
n <- nrow(subjData)


######################################
### JLF cortical thickness regions ###
######################################

#Get variable names
ctRegions <- names(subjData)[grep("mprage_jlf_ct",names(subjData))]

##Run models (all three groups included).
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
ctGam <- lapply(ctRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + white + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
ctGamSumm <- lapply(ctGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
ctAnova <- lapply(ctGam, anova)

#Pull the p-values
ct_p <- sapply(ctAnova, function(v) v$pTerms.pv[4])

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
ct_omnibus <- cbind(ctNames,ct_p_fdr_round)

#Trim table to only regions that passed FDR correction
ct_omnibus_signif <- ct_omnibus[ct_p_fdr<0.05,]

#Save omnibus results as a .csv file
write.table(ct_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_ct_race.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs. S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
ctGam_reordered <- lapply(ctRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
ctGamSumm_reordered <- lapply(ctGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
ct_S1vsTd <- sapply(ctGam, function(v) summary(v)$p.table[5,4])
ct_S2vsTd <- sapply(ctGam, function(v) summary(v)$p.table[6,4])
ct_S1vsS2 <- sapply(ctGam_reordered, function(v) summary(v)$p.table[5,4])

#Combine the pairwise p values
ct_pairs <- cbind(ct_S1vsTd,ct_S2vsTd,ct_S1vsS2)

#Convert to data frame
ct_pairs <- as.data.frame(ct_pairs)

#Add row names for ease of viewing
rownames(ct_pairs) <- ctRegions

#Print original p-values to three decimal places
ct_pairs_round <- round(ct_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
ct_fdrTable<-as.data.frame(matrix(nrow=98,ncol=3))
colnames(ct_fdrTable)[1]<-"ct_S1vsTd_pfdr"
colnames(ct_fdrTable)[2]<-"ct_S2vsTd_pfdr"
colnames(ct_fdrTable)[3]<-"ct_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(ct_pairs)) {
    row <- ct_pairs[i,]
    ct_fdrTable[i,] <- p.adjust(ct_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
ct_fdrTable_round <- round(ct_fdrTable,3)

#Add region names
ct_pairwise <- cbind(ctNames,ct_fdrTable_round)

#Trim table to only regions that passed FDR correction for the omnibus test
ct_pairwise_signif <- ct_pairwise[ct_p_fdr<0.05,]

#Save the pairwise results as a .csv file
write.table(ct_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_ct_race.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##################
### Average CT ###
##################

#Create a mean ct variable that includes all significant ct regions
subjData$averageCT_signifRegions <- rowMeans(subset(subjData, select = ctRegions))

#Run model
avgCtGam <- gam(averageCT_signifRegions ~ s(age) + sex + averageManualRating + white + Hydra_k2, method="REML", data=subjData)

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
avgCtGam2 <- gam(averageCT_signifRegions ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
avgCtAnova <- anova(avgCtGam)

##Pairwise comparisons
#Pull uncorrected p-values
avgCt_S1vsTd <- summary(avgCtGam)$p.table[5,4]
avgCt_S2vsTd <- summary(avgCtGam)$p.table[6,4]
avgCt_S1vsS2 <- summary(avgCtGam2)$p.table[5,4]

#Combine the pairwise p values
avgCt_pairs <- cbind(avgCt_S1vsTd,avgCt_S2vsTd,avgCt_S1vsS2)

#Convert to data frame
avgCt_pairs <- as.data.frame(avgCt_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
avgCt_fdr <- p.adjust(avgCt_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
avgCt_fdr_round <- round(avgCt_fdr,3)


##########################
### JLF volume regions ###
##########################

#Get variable names
volRegions_all <- subjData[c(grep("mprage_jlf_vol",names(subjData)))]
volRegions_short <- volRegions_all[,-grep("Vent|Brain_Stem|Cerebell|Cerebral_White_Matter|CSF|Lobe_WM",names(volRegions_all))]
volRegions <- names(volRegions_short)

##Run models (all three groups included).
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
volGam <- lapply(volRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + white + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
volGamSumm <- lapply(volGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
volAnova <- lapply(volGam, anova)

#Pull the p-values
vol_p <- sapply(volAnova, function(v) v$pTerms.pv[4])

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
vol_omnibus <- cbind(volNames,vol_p_fdr_round)

#Trim table to only regions that passed FDR correction
vol_omnibus_signif <- vol_omnibus[vol_p_fdr<0.05,]

#Save omnibus results as a .csv file
write.table(vol_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_vol_race.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
volGam_reordered <- lapply(volRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
volGamSumm_reordered <- lapply(volGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
vol_S1vsTd <- sapply(volGam, function(v) summary(v)$p.table[5,4])
vol_S2vsTd <- sapply(volGam, function(v) summary(v)$p.table[6,4])
vol_S1vsS2 <- sapply(volGam_reordered, function(v) summary(v)$p.table[5,4])

#Combine the pairwise p values
vol_pairs <- cbind(vol_S1vsTd,vol_S2vsTd,vol_S1vsS2)

#Convert to data frame
vol_pairs <- as.data.frame(vol_pairs)

#Add row names for ease of viewing
rownames(vol_pairs) <- volRegions

#Print original p-values to three decimal places
vol_pairs_round <- round(vol_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
vol_fdrTable<-as.data.frame(matrix(nrow=112,ncol=3))
colnames(vol_fdrTable)[1]<-"vol_S1vsTd_pfdr"
colnames(vol_fdrTable)[2]<-"vol_S2vsTd_pfdr"
colnames(vol_fdrTable)[3]<-"vol_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(vol_pairs)) {
    row <- vol_pairs[i,]
    vol_fdrTable[i,] <- p.adjust(vol_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
vol_fdrTable_round <- round(vol_fdrTable,3)

#Add region names
vol_pairwise <- cbind(volNames,vol_fdrTable_round)

#Trim table to only regions that passed FDR correction for the omnibus test
vol_pairwise_signif <- vol_pairwise[vol_p_fdr<0.05,]

#Save the pairwise results as a .csv file
write.table(vol_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_vol_race.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


################################
### Total gray matter volume ###
################################

#Divide total gray matter volume by 1000 to change the units from cubic millimeters (mm3) to cubic centimeters (cc3); 1 cc3 = 1,000 mm3
subjData$mprage_antsCT_vol_GrayMatter <- subjData$mprage_antsCT_vol_GrayMatter/1000

#Run model
totGrayGam <- gam(mprage_antsCT_vol_GrayMatter ~ s(age) + sex + averageManualRating + white + Hydra_k2, method="REML", data=subjData)

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
totGrayGam2 <- gam(mprage_antsCT_vol_GrayMatter ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
totGrayAnova <- anova(totGrayGam)

##Pairwise comparisons
#Pull uncorrected p-values
totGray_S1vsTd <- summary(totGrayGam)$p.table[5,4]
totGray_S2vsTd <- summary(totGrayGam)$p.table[6,4]
totGray_S1vsS2 <- summary(totGrayGam2)$p.table[5,4]

#Combine the pairwise p values
totGray_pairs <- cbind(totGray_S1vsTd,totGray_S2vsTd,totGray_S1vsS2)

#Convert to data frame
totGray_pairs <- as.data.frame(totGray_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
totGray_fdr <- p.adjust(totGray_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
totGray_fdr_round <- round(totGray_fdr,3)


###########
### ICV ###
###########

#Divide ICV by 1000 to change the units from cubic millimeters (mm3) to cubic centimeters (cc3); 1 cc3 = 1,000 mm3
subjData$ICV <- subjData$ICV/1000

#Run model
icvGam <- gam(ICV ~ s(age) + sex + averageManualRating + white + Hydra_k2, method="REML", data=subjData)

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
icvGam2 <- gam(ICV ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
icvAnova <- anova(icvGam)

##Pairwise comparisons
#Pull uncorrected p-values
icv_S1vsTd <- summary(icvGam)$p.table[5,4]
icv_S2vsTd <- summary(icvGam)$p.table[6,4]
icv_S1vsS2 <- summary(icvGam2)$p.table[5,4]

#Combine the pairwise p values
icv_pairs <- cbind(icv_S1vsTd,icv_S2vsTd,icv_S1vsS2)

#Convert to data frame
icv_pairs <- as.data.frame(icv_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
icv_fdr <- p.adjust(icv_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
icv_fdr_round <- round(icv_fdr,3)


###########
### TBV ###
###########

#Divide TBV by 1000 to change the units from cubic millimeters (mm3) to cubic centimeters (cc3); 1 cc3 = 1,000 mm3
subjData$mprage_antsCT_vol_TBV <- subjData$mprage_antsCT_vol_TBV/1000

#Run model
tbvGam <- gam(mprage_antsCT_vol_TBV ~ s(age) + sex + averageManualRating + white + Hydra_k2, method="REML", data=subjData)

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
tbvGam2 <- gam(mprage_antsCT_vol_TBV ~ s(age) + sex + averageManualRating + white + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
tbvAnova <- anova(tbvGam)

##Pairwise comparisons
#Pull uncorrected p-values
tbv_S1vsTd <- summary(tbvGam)$p.table[5,4]
tbv_S2vsTd <- summary(tbvGam)$p.table[6,4]
tbv_S1vsS2 <- summary(tbvGam2)$p.table[5,4]

#Combine the pairwise p values
tbv_pairs <- cbind(tbv_S1vsTd,tbv_S2vsTd,tbv_S1vsS2)

#Convert to data frame
tbv_pairs <- as.data.frame(tbv_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
tbv_fdr <- p.adjust(tbv_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
tbv_fdr_round <- round(tbv_fdr,3)
