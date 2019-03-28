##################################################
#### GROUP ANALYSES WITH VOLUME/CT JLF SAMPLE ####
##################################################

###############################
### Load data and libraries ###
###############################

subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_subjData.rds")

#Load libraries ("plyr" needs to be loaded before "dplyr")
library(plyr)
library(mgcv)
library(dplyr)
library(varhandle)
library(psych)
library(schoolmath)

#Total sample size
n <- nrow(subjData)


########################
### Rename variables ###
########################

#The right accumbens is named "R_Accumbens_Area" in the dataset but "R_Accumbens" in the JLF labels csv. Rename to match.
names(subjData)[names(subjData) == 'mprage_jlf_ct_R_Accumbens_Area'] <- 'mprage_jlf_ct_R_Accumbens'
names(subjData)[names(subjData) == 'mprage_jlf_ct_L_Accumbens_Area'] <- 'mprage_jlf_ct_L_Accumbens'

names(subjData)[names(subjData) == 'mprage_jlf_vol_R_Accumbens_Area'] <- 'mprage_jlf_vol_R_Accumbens'
names(subjData)[names(subjData) == 'mprage_jlf_vol_L_Accumbens_Area'] <- 'mprage_jlf_vol_L_Accumbens'


####################
### Descriptives ###
####################

#How many subjects in each group (0=TD, 1=subtype1, 2=subtype2)
SampleSize_table <- table(subjData$Hydra_k2) 

#Group means (0=TD, 1=subtype1, 2=subtype2)
meanSdAge <- ddply(subjData,~Hydra_k2,summarise,mean=mean(age),sd=sd(age))
meanSdMedu <- ddply(subjData,~Hydra_k2,summarise,mean=mean(medu1,na.rm=TRUE),sd=sd(medu1,na.rm=TRUE))

#Percentage of females (0=male, 1=female) (sex can't be a factor and needs to be 0s and 1s)
subjData_temp <- subjData
subjData_temp$sex <- unfactor(subjData_temp$sex)
percentFemale <- ddply(subjData_temp,~Hydra_k2,summarise,mean=mean(sex))

#Or to see the number of males and females in each group
sexTable <- table(Hydra_k2 = subjData$Hydra_k2, Sex = subjData$sex)


#############################################
#### Lifetime prevalence psychopathology ####
#############################################

#Typically Developing (N and percent)
Td_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Td, na.rm=TRUE))
Td_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Td, na.rm=TRUE)/1141)

#ADHD Diagnosis
Add_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Add, na.rm=TRUE))
Add_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Add, na.rm=TRUE)/1141)

#Agoraphobia Diagnosis
Agr_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Agr, na.rm=TRUE))
Agr_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Agr, na.rm=TRUE)/1141)

#Anorexia Diagnosis
Ano_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Ano, na.rm=TRUE))
Ano_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Ano, na.rm=TRUE)/1141)

#Bulimia Diagnosis
Bul_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Bul, na.rm=TRUE))
Bul_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Bul, na.rm=TRUE)/1141)

#Conduct Disorder Diagnosis
Con_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Con, na.rm=TRUE))
Con_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Con, na.rm=TRUE)/1141)

#Generalized Anxiety Disorder Diagnosis
Gad_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Gad, na.rm=TRUE))
Gad_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Gad, na.rm=TRUE)/1141)

#Major Depression Diagnosis
Mdd_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Mdd, na.rm=TRUE))
Mdd_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Mdd, na.rm=TRUE)/1141)

#Mania Diagnosis
Man_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Man, na.rm=TRUE))
Man_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Man, na.rm=TRUE)/1141)

#OCD Diagnosis
Ocd_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Ocd, na.rm=TRUE))
Ocd_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Ocd, na.rm=TRUE)/1141)

#ODD Diagnosis
Odd_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Odd, na.rm=TRUE))
Odd_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Odd, na.rm=TRUE)/1141)

#Panic Diagnosis
Pan_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Pan, na.rm=TRUE))
Pan_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Pan, na.rm=TRUE)/1141)

#Psychosis spectrum Diagnosis
Ps_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Ps, na.rm=TRUE))
Ps_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Ps, na.rm=TRUE)/1141)

#PTSD Diagnosis
Ptd_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Ptd, na.rm=TRUE))
Ptd_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Ptd, na.rm=TRUE)/1141)

#Seperation Anxiety Diagnosis
Sep_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Sep, na.rm=TRUE))
Sep_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Sep, na.rm=TRUE)/1141)

#Social Phobia Diagnosis
Soc_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Soc, na.rm=TRUE))
Soc_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Soc, na.rm=TRUE)/1141)

#Specific Phobia Diagnosis
Sph_n <- ddply(subjData,~Hydra_k2,summarise,N=sum(Sph, na.rm=TRUE))
Sph_percent <- ddply(subjData,~Hydra_k2,summarise,percent=sum(Sph, na.rm=TRUE)/1141)


##############
### Gender ###
##############

##Run model
genderLm <- glm(sex ~ Hydra_k2, data=subjData, family=binomial)

#Look at model summary
genderLmSumm <- summary(genderLm)


##Pass the results to the anova() function to fit a Chi-square test
#This is your omnibus Chi-square test (tells you if sex differs by group)
#df1=k-1, df2=n-k
genderAnova <- anova(genderLm, test="Chisq")


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
genderLm_reordered <- glm(sex ~ Hydra_k2_reordered, data=subjData, family=binomial)

#Look at model summary
genderLmSumm_reordered <- summary(genderLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
gender_S1vsTd <- summary(genderLm)$coefficients[2,4]
gender_S2vsTd <- summary(genderLm)$coefficients[3,4]
gender_S1vsS2 <- summary(genderLm_reordered)$coefficients[2,4]

#Combine the pairwise p values
gender_pairs <- cbind(gender_S1vsTd,gender_S2vsTd,gender_S1vsS2)

#Convert to data frame
gender_pairs <- as.data.frame(gender_pairs)

#Print original p-values to three decimal places
gender_pairs_round <- round(gender_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
gender_fdr <- p.adjust(gender_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
gender_fdr_round <- round(gender_fdr,3)


###########
### Age ###
###########

##Run model
ageLm <- lm(age ~ Hydra_k2, data=subjData)

#Look at model summary
ageLmSumm <- summary(ageLm)


##Pass the results to the anova() function to fit ANOVAs
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
ageAnova <- anova(ageLm)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
ageLm_reordered <- lm(age ~ Hydra_k2_reordered, data=subjData)

#Look at model summary
ageLmSumm_reordered <- summary(ageLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
age_S1vsTd <- summary(ageLm)$coefficients[2,4]
age_S2vsTd <- summary(ageLm)$coefficients[3,4]
age_S1vsS2 <- summary(ageLm_reordered)$coefficients[2,4]

#Combine the pairwise p values
age_pairs <- cbind(age_S1vsTd,age_S2vsTd,age_S1vsS2)

#Convert to data frame
age_pairs <- as.data.frame(age_pairs)

#Print original p-values to three decimal places
age_pairs_round <- round(age_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
age_fdr <- p.adjust(age_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
age_fdr_round <- round(age_fdr,3)


##Follow-up on significant results
#Means and SDs
ageDescriptives <- tapply(subjData$age, subjData$Hydra_k2, describe)


###################################
### Maternal level of education ###
###################################

##Run model
meduLm <- lm(medu1 ~ Hydra_k2, data=subjData)

#Look at model summary
meduLmSumm <- summary(meduLm)


##Pass the results to the anova() function to fit ANOVAs
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
meduAnova <- anova(meduLm)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
meduLm_reordered <- lm(medu1 ~ Hydra_k2_reordered, data=subjData)

#Look at model summary
meduLmSumm_reordered <- summary(meduLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
medu_S1vsTd <- summary(meduLm)$coefficients[2,4]
medu_S2vsTd <- summary(meduLm)$coefficients[3,4]
medu_S1vsS2 <- summary(meduLm_reordered)$coefficients[2,4]

#Combine the pairwise p values
medu_pairs <- cbind(medu_S1vsTd,medu_S2vsTd,medu_S1vsS2)

#Convert to data frame
medu_pairs <- as.data.frame(medu_pairs)

#Print original p-values to three decimal places
medu_pairs_round <- round(medu_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
medu_fdr <- p.adjust(medu_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
medu_fdr_round <- round(medu_fdr,3)


##Follow-up on significant results
#Means and SDs
meduDescriptives <- tapply(subjData$medu1, subjData$Hydra_k2, describe)


#################################
### Psychopathology bifactors ###
#################################

#Get bifactor variable names
bifactors <- names(subjData)[grep("4factor",names(subjData))]

#Run models (all three groups included)
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
biLm <- lapply(bifactors, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
biLmSumm <- lapply(biLm, summary)


##Pass the results to the anova() function to fit ANOVAs
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
biAnova <- lapply(biLm, anova)

#Pull the p-values
bi_p <- sapply(biAnova, function(v) v$"Pr(>F)"[3])

#Convert to data frame
bi_p <- as.data.frame(bi_p)

#Add row names for ease of viewing
rownames(bi_p) <- bifactors

#Print original p-values to three decimal places
bi_p_round <- round(bi_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values
bi_p_fdr <- p.adjust(bi_p[,1],method="fdr")

#Convert to data frame
bi_p_fdr <- as.data.frame(bi_p_fdr)

#Print fdr-corrected p-values to three decimal places
bi_p_fdr_round <- round(bi_p_fdr,3)

#Add names
biNames <- as.data.frame(bifactors)
bi_omnibus <- cbind(biNames,bi_p_fdr_round)

#Trim table to only variables that passed FDR correction
bi_omnibus_signif <- bi_omnibus[bi_p_fdr<0.05,]


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
biLm_reordered <- lapply(bifactors, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2_reordered, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
biLmSumm_reordered <- lapply(biLm_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
bi_S1vsTd <- sapply(biLm, function(v) summary(v)$coefficients[4,4])
bi_S2vsTd <- sapply(biLm, function(v) summary(v)$coefficients[5,4])
bi_S1vsS2 <- sapply(biLm_reordered, function(v) summary(v)$coefficients[4,4])

#Combine the pairwise p values
bi_pairs <- cbind(bi_S1vsTd,bi_S2vsTd,bi_S1vsS2)

#Convert to data frame
bi_pairs <- as.data.frame(bi_pairs)

#Add row names for ease of viewing
rownames(bi_pairs) <- bifactors

#Print original p-values to three decimal places
bi_pairs_round <- round(bi_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
bi_fdrTable<-as.data.frame(matrix(nrow=5,ncol=3))
colnames(bi_fdrTable)[1]<-"bi_S1vsTd_pfdr"
colnames(bi_fdrTable)[2]<-"bi_S2vsTd_pfdr"
colnames(bi_fdrTable)[3]<-"bi_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(bi_pairs)) {
    row <- bi_pairs[i,]
    bi_fdrTable[i,] <- p.adjust(bi_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
bi_fdrTable_round <- round(bi_fdrTable,3)

#Add names
bi_pairwise <- cbind(biNames,bi_fdrTable_round)


##Follow-up on significant results
#Means and SDs
overallDescriptives <- tapply(subjData$overall_psychopathology_4factorv2, subjData$Hydra_k2, describe)
extDescriptives <- tapply(subjData$externalizing_4factorv2, subjData$Hydra_k2, describe)
fearDescriptives <- tapply(subjData$phobias_4factorv2, subjData$Hydra_k2, describe)


########################
### Overall Accuracy ###
########################

##Run model
cogLm <- lm(Overall_Accuracy ~ age + sex + Hydra_k2, data=subjData)

#Look at model summary
cogLmSumm <- summary(cogLm)


##Pass the results to the anova() function to fit ANOVAs
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
cogAnova <- anova(cogLm)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
cogLm_reordered <- lm(Overall_Accuracy ~ age + sex + Hydra_k2_reordered, data=subjData)

#Look at model summary
cogLmSumm_reordered <- summary(cogLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
cog_S1vsTd <- summary(cogLm)$coefficients[4,4]
cog_S2vsTd <- summary(cogLm)$coefficients[5,4]
cog_S1vsS2 <- summary(cogLm_reordered)$coefficients[4,4]

#Combine the pairwise p values
cog_pairs <- cbind(cog_S1vsTd,cog_S2vsTd,cog_S1vsS2)

#Convert to data frame
cog_pairs <- as.data.frame(cog_pairs)

#Print original p-values to three decimal places
cog_pairs_round <- round(cog_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
cog_fdr <- p.adjust(cog_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
cog_fdr_round <- round(cog_fdr,3)


##Follow-up on significant results
#Means and SDs
overallAccDescriptives <- tapply(subjData$Overall_Accuracy, subjData$Hydra_k2, describe)


##################################
### Cognitive Accuracy Factors ###
##################################

#Get correlated traits variable names
AccFactors <- names(subjData)[grep("^F1_Exec_Comp_Res_Accuracy$|^F2_Social_Cog_Accuracy$|^F3_Memory_Accuracy$",names(subjData))]

#Run models (all three groups included)
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
AccFacLm <- lapply(AccFactors, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
AccFacLmSumm <- lapply(AccFacLm, summary)


##Pass the results to the anova() function to fit ANOVAs
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
AccFacAnova <- lapply(AccFacLm, anova)

#Pull the p-values
AccFac_p <- sapply(AccFacAnova, function(v) v$"Pr(>F)"[3])

#Convert to data frame
AccFac_p <- as.data.frame(AccFac_p)

#Add row names for ease of viewing
rownames(AccFac_p) <- AccFactors

#Print original p-values to three decimal places
AccFac_p_round <- round(AccFac_p,3)


##FDR correction across all omnibus anova tests
#FDR correct p-values
AccFac_p_fdr <- p.adjust(AccFac_p[,1],method="fdr")

#Convert to data frame
AccFac_p_fdr <- as.data.frame(AccFac_p_fdr)

#Print fdr-corrected p-values to three decimal places
AccFac_p_fdr_round <- round(AccFac_p_fdr,3)

#Add names
AccFacNames <- as.data.frame(AccFactors)
AccFac_omnibus <- cbind(AccFacNames,AccFac_p_fdr_round)

#Trim table to only variables that passed FDR correction
AccFac_omnibus_signif <- AccFac_omnibus[AccFac_p_fdr<0.05,]


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
AccFacLm_reordered <- lapply(AccFactors, function(x) {
  lm(substitute(i ~ age + sex + Hydra_k2_reordered, list(i = as.name(x))), data=subjData)
})

#Look at model summaries
AccFacLmSumm_reordered <- lapply(AccFacLm_reordered, summary)

##Pairwise comparisons
#Pull uncorrected p-values
AccFac_S1vsTd <- sapply(AccFacLm, function(v) summary(v)$coefficients[4,4])
AccFac_S2vsTd <- sapply(AccFacLm, function(v) summary(v)$coefficients[5,4])
AccFac_S1vsS2 <- sapply(AccFacLm_reordered, function(v) summary(v)$coefficients[4,4])

#Combine the pairwise p values
AccFac_pairs <- cbind(AccFac_S1vsTd,AccFac_S2vsTd,AccFac_S1vsS2)

#Convert to data frame
AccFac_pairs <- as.data.frame(AccFac_pairs)

#Add row names for ease of viewing
rownames(AccFac_pairs) <- AccFactors

#Print original p-values to three decimal places
AccFac_pairs_round <- round(AccFac_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
#Create an empty table for the fdr results
AccFac_fdrTable<-as.data.frame(matrix(nrow=3,ncol=3))
colnames(AccFac_fdrTable)[1]<-"AccFac_S1vsTd_pfdr"
colnames(AccFac_fdrTable)[2]<-"AccFac_S2vsTd_pfdr"
colnames(AccFac_fdrTable)[3]<-"AccFac_S1vsS2_pfdr"

#FDR correct across rows
for(i in 1:nrow(AccFac_pairs)) {
    row <- AccFac_pairs[i,]
    AccFac_fdrTable[i,] <- p.adjust(AccFac_pairs[i,],method="fdr")
}

#Print fdr-corrected p-values to three decimal places
AccFac_fdrTable_round <- round(AccFac_fdrTable,3)

#Add names
AccFac_pairwise <- cbind(AccFacNames,AccFac_fdrTable_round)


##Follow-up on significant results
#Means and SDs
Comp_Res_Acc_Descriptives <- tapply(subjData$F1_Exec_Comp_Res_Accuracy, subjData$Hydra_k2, describe)
Social_Cog_Acc_Descriptives <- tapply(subjData$F2_Social_Cog_Accuracy, subjData$Hydra_k2, describe)
Memory_Acc_Descriptives <- tapply(subjData$F3_Memory_Accuracy, subjData$Hydra_k2, describe)


#####################
#### WRAT scores ####
#####################

##Run model
wratLm <- lm(wrat4CrStd ~ age + sex + Hydra_k2, data=subjData)

#Look at model summary
wratLmSumm <- summary(wratLm)


##Pass the results to the anova() function to fit ANOVAs
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
wratAnova <- anova(wratLm)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
wratLm_reordered <- lm(wrat4CrStd ~ age + sex + Hydra_k2_reordered, data=subjData)

#Look at model summary
wratLmSumm_reordered <- summary(wratLm_reordered)


##Pairwise comparisons
#Pull uncorrected p-values
wrat_S1vsTd <- summary(wratLm)$coefficients[4,4]
wrat_S2vsTd <- summary(wratLm)$coefficients[5,4]
wrat_S1vsS2 <- summary(wratLm_reordered)$coefficients[4,4]

#Combine the pairwise p values
wrat_pairs <- cbind(wrat_S1vsTd,wrat_S2vsTd,wrat_S1vsS2)

#Convert to data frame
wrat_pairs <- as.data.frame(wrat_pairs)

#Print original p-values to three decimal places
wrat_pairs_round <- round(wrat_pairs,3)


##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
wrat_fdr <- p.adjust(wrat_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
wrat_fdr_round <- round(wrat_fdr,3)


##Follow-up on significant results
#Means and SDs
wratDescriptives <- tapply(subjData$wrat4CrStd, subjData$Hydra_k2, describe)


##############################
### JLF cortical thickness ###
##############################

#Get variable names
ctRegions <- names(subjData)[grep("mprage_jlf_ct",names(subjData))]

##Run models (all three groups included).
#TD is the comparison group in Hydra_k2. This model will provide estimates and p-values for TD vs. S1 and TD vs. S2.
ctGam <- lapply(ctRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
ctGamSumm <- lapply(ctGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
ctAnova <- lapply(ctGam, anova)

#Pull the p-values
ct_p <- sapply(ctAnova, function(v) v$pTerms.pv[3])

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
write.table(ct_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_ct.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs. S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
ctGam_reordered <- lapply(ctRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
ctGamSumm_reordered <- lapply(ctGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
ct_S1vsTd <- sapply(ctGam, function(v) summary(v)$p.table[4,4])
ct_S2vsTd <- sapply(ctGam, function(v) summary(v)$p.table[5,4])
ct_S1vsS2 <- sapply(ctGam_reordered, function(v) summary(v)$p.table[4,4])

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
write.table(ct_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_ct.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Significant S1 vs. S2 t-values
#Pull t-values
ct_t_S1vsS2 <- sapply(ctGam_reordered, function(v) summary(v)$p.table[4,3])

#Convert to data frame
ct_t_S1vsS2 <- as.data.frame(ct_t_S1vsS2)

#Print t-values to two decimal places
ct_t_S1vsS2_round <- round(ct_t_S1vsS2,2)

#Add region names
ct_t_S1vsS2_names <- cbind(ctNames,ct_t_S1vsS2_round)

#Trim table to only regions that passed FDR correction for the omnibus test
ct_t_S1vsS2_omnibus <- ct_t_S1vsS2_names[ct_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs S2 pairwise test
ct_t_S1vsS2_signif <- ct_t_S1vsS2_omnibus[ct_pairwise_signif$ct_S1vsS2_pfdr<0.05,]

#Direction of results
ct_t_S1vsS2_neg <- is.negative(ct_t_S1vsS2_signif$ct_t_S1vsS2)


##Create table of region names and t-values for mapping in fslview
#Load JLF labels and corresponding index numbers
#These are taken from https://github.com/PennBBL/jlfVisualizer/blob/master/data/jlf_lookupWithWM.csv
JLF.labels <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/images/jlf_lookupWithWM.csv", header=TRUE)

#Rename variables for merging
ct_t_S1vsS2_renamed <- rename(ct_t_S1vsS2_signif, ROI_NAME = ctRegions)

#Remove the "mprage_jlf_ct_" prefix to match with the JLF labels
ct_t_S1vsS2_renamed[] <- lapply(ct_t_S1vsS2_renamed, function(x) gsub("mprage_jlf_ct_", "", x))

#Merge to add index numbers that correspond to the significant regions
ct_merged_S1vsS2 <-merge(ct_t_S1vsS2_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
ct_merged_S1vsS2$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only	makes sense if all your results are in the same direction (all negative).
ct_merged_S1vsS2$ct_t_S1vsS2 <- as.numeric(ct_merged_S1vsS2$ct_t_S1vsS2)
ct_merged_S1vsS2$ct_t_S1vsS2 <- abs(ct_merged_S1vsS2$ct_t_S1vsS2)

#Save as a .csv
write.table(ct_merged_S1vsS2, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_ct_S1vsS2.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Significant S1 vs. TD t-values
#Pull t-values
ct_t_S1vsTd <- sapply(ctGam, function(v) summary(v)$p.table[4,3])

#Convert to data frame
ct_t_S1vsTd <- as.data.frame(ct_t_S1vsTd)

#Print t-values to two decimal places
ct_t_S1vsTd_round <- round(ct_t_S1vsTd,2)

#Add region names
ct_t_S1vsTd_names <- cbind(ctNames,ct_t_S1vsTd_round)

#Trim table to only regions that passed FDR correction for the omnibus test
ct_t_S1vsTd_omnibus <- ct_t_S1vsTd_names[ct_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs TD pairwise test
ct_t_S1vsTd_signif <- ct_t_S1vsTd_omnibus[ct_pairwise_signif$ct_S1vsTd_pfdr<0.05,]

#Direction of results
ct_t_S1vsTd_neg <- is.negative(ct_t_S1vsTd_signif$ct_t_S1vsTd)


##Create table of region names and t-values for mapping in fslview
#Rename variables for merging
ct_t_S1vsTd_renamed <- rename(ct_t_S1vsTd_signif, ROI_NAME = ctRegions)

#Remove the "mprage_jlf_ct_" prefix to match with the JLF labels
ct_t_S1vsTd_renamed[] <- lapply(ct_t_S1vsTd_renamed, function(x) gsub("mprage_jlf_ct_", "", x))

#Merge to add index numbers that correspond to the significant regions
ct_merged_S1vsTd <-merge(ct_t_S1vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
ct_merged_S1vsTd$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
ct_merged_S1vsTd$ct_t_S1vsTd <- as.numeric(ct_merged_S1vsTd$ct_t_S1vsTd)
ct_merged_S1vsTd$ct_t_S1vsTd <- abs(ct_merged_S1vsTd$ct_t_S1vsTd)

#Save as a .csv
write.table(ct_merged_S1vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_ct_S1vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Significant S2 vs. TD t-values
#Pull t-values
ct_t_S2vsTd <- sapply(ctGam, function(v) summary(v)$p.table[5,3])

#Convert to data frame
ct_t_S2vsTd <- as.data.frame(ct_t_S2vsTd)

#Print t-values to two decimal places
ct_t_S2vsTd_round <- round(ct_t_S2vsTd,2)

#Add region names
ct_t_S2vsTd_names <- cbind(ctNames,ct_t_S2vsTd_round)

#Trim table to only regions that passed FDR correction for the omnibus test
ct_t_S2vsTd_omnibus <- ct_t_S2vsTd_names[ct_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S2 vs Td pairwise test
ct_t_S2vsTd_signif <- ct_t_S2vsTd_omnibus[ct_pairwise_signif$ct_S2vsTd_pfdr<0.05,]

#Direction of results
ct_t_S2vsTd_neg <- is.negative(ct_t_S2vsTd_signif$ct_t_S2vsTd)


##Create table of region names and t-values for mapping in fslview
#Rename variables for merging
ct_t_S2vsTd_renamed <- rename(ct_t_S2vsTd_signif, ROI_NAME = ctRegions)

#Remove the "mprage_jlf_ct_" prefix to match with the JLF labels
ct_t_S2vsTd_renamed[] <- lapply(ct_t_S2vsTd_renamed, function(x) gsub("mprage_jlf_ct_", "", x))

#Merge to add index numbers that correspond to the significant regions
ct_merged_S2vsTd <-merge(ct_t_S2vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
ct_merged_S2vsTd$ROI_NAME <- NULL

#Save as a .csv
write.table(ct_merged_S2vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_ct_S2vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


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
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
volGamSumm <- lapply(volGam, summary)


##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#These are your omnibus ANOVA tests (tells you if the three groups are significantly different)
volAnova <- lapply(volGam, anova)

#Pull the p-values
vol_p <- sapply(volAnova, function(v) v$pTerms.pv[3])

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
write.table(vol_omnibus_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_omnibus_vol.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
#This model will provide estimates and p-values for S1-S2.
volGam_reordered <- lapply(volRegions, function(x) {
  gam(substitute(i ~ s(age) + sex + averageManualRating + Hydra_k2_reordered, list(i = as.name(x))), method="REML", data=subjData)
})

#Look at model summaries
volGamSumm_reordered <- lapply(volGam_reordered, summary)


##Pairwise comparisons
#Pull uncorrected p-values
vol_S1vsTd <- sapply(volGam, function(v) summary(v)$p.table[4,4])
vol_S2vsTd <- sapply(volGam, function(v) summary(v)$p.table[5,4])
vol_S1vsS2 <- sapply(volGam_reordered, function(v) summary(v)$p.table[4,4])

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
write.table(vol_pairwise_signif, file = "/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/pfdr_pairwise_vol.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


##Significant S1 vs. S2 t-values
#Pull t-values
vol_t_S1vsS2 <- sapply(volGam_reordered, function(v) summary(v)$p.table[4,3])

#Convert to data frame
vol_t_S1vsS2 <- as.data.frame(vol_t_S1vsS2)

#Print t-values to two decimal places
vol_t_S1vsS2_round <- round(vol_t_S1vsS2,2)

#Add region names
vol_t_S1vsS2_names <- cbind(volNames,vol_t_S1vsS2_round)

#Trim table to only regions that passed FDR correction for the omnibus test
vol_t_S1vsS2_omnibus <- vol_t_S1vsS2_names[vol_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs S2 pairwise test
vol_t_S1vsS2_signif <- vol_t_S1vsS2_omnibus[vol_pairwise_signif$vol_S1vsS2_pfdr<0.05,]

#Direction of results
vol_t_S1vsS2_neg <- is.negative(vol_t_S1vsS2_signif$vol_t_S1vsS2)


##Create table of region names and t-values for mapping in fslview
#Rename variables for merging
vol_t_S1vsS2_renamed <- rename(vol_t_S1vsS2_signif, ROI_NAME = volRegions)

#Remove the "mprage_jlf_vol_" prefix to match with the JLF labels
vol_t_S1vsS2_renamed[] <- lapply(vol_t_S1vsS2_renamed, function(x) gsub("mprage_jlf_vol_", "", x))

#Merge to add index numbers that correspond to the significant regions
vol_merged_S1vsS2 <-merge(vol_t_S1vsS2_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
vol_merged_S1vsS2$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
vol_merged_S1vsS2$vol_t_S1vsS2 <- as.numeric(vol_merged_S1vsS2$vol_t_S1vsS2)
vol_merged_S1vsS2$vol_t_S1vsS2 <- abs(vol_merged_S1vsS2$vol_t_S1vsS2)

#Save as a .csv
write.table(vol_merged_S1vsS2, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_vol_S1vsS2.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")


##Significant S1 vs. TD t-values
#Pull t-values
vol_t_S1vsTd <- sapply(volGam, function(v) summary(v)$p.table[4,3])

#Convert to data frame
vol_t_S1vsTd <- as.data.frame(vol_t_S1vsTd)

#Print t-values to two decimal places
vol_t_S1vsTd_round <- round(vol_t_S1vsTd,2)

#Add region names
vol_t_S1vsTd_names <- cbind(volNames,vol_t_S1vsTd_round)

#Trim table to only regions that passed FDR correction for the omnibus test
vol_t_S1vsTd_omnibus <- vol_t_S1vsTd_names[vol_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S1 vs TD pairwise test
vol_t_S1vsTd_signif <- vol_t_S1vsTd_omnibus[vol_pairwise_signif$vol_S1vsTd_pfdr<0.05,]

#Direction of results
vol_t_S1vsTd_neg <- is.negative(vol_t_S1vsTd_signif$vol_t_S1vsTd)


##Create table of region names and t-values for mapping in fslview
#Rename variables for merging
vol_t_S1vsTd_renamed <- rename(vol_t_S1vsTd_signif, ROI_NAME = volRegions)

#Remove the "mprage_jlf_vol_" prefix to match with the JLF labels
vol_t_S1vsTd_renamed[] <- lapply(vol_t_S1vsTd_renamed, function(x) gsub("mprage_jlf_vol_", "", x))

#Merge to add index numbers that correspond to the significant regions
vol_merged_S1vsTd <-merge(vol_t_S1vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
vol_merged_S1vsTd$ROI_NAME <- NULL

#Make all the values positive instead of negative just so we can visualize the strongest regions in light blue instead of dark blue in fslview
#Note this only makes sense if all your results are in the same direction (all negative).
vol_merged_S1vsTd$vol_t_S1vsTd <- as.numeric(vol_merged_S1vsTd$vol_t_S1vsTd)
vol_merged_S1vsTd$vol_t_S1vsTd <- abs(vol_merged_S1vsTd$vol_t_S1vsTd)

#Save as a .csv
write.table(vol_merged_S1vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_vol_S1vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")



##Significant S2 vs. TD t-values
#Pull t-values
vol_t_S2vsTd <- sapply(volGam, function(v) summary(v)$p.table[5,3])

#Convert to data frame
vol_t_S2vsTd <- as.data.frame(vol_t_S2vsTd)

#Print t-values to two decimal places
vol_t_S2vsTd_round <- round(vol_t_S2vsTd,2)

#Add region names
vol_t_S2vsTd_names <- cbind(volNames,vol_t_S2vsTd_round)

#Trim table to only regions that passed FDR correction for the omnibus test
vol_t_S2vsTd_omnibus <- vol_t_S2vsTd_names[vol_p_fdr<0.05,]

#Trim table to only regions that passed FDR correction for the S2 vs TD pairwise test
vol_t_S2vsTd_signif <- vol_t_S2vsTd_omnibus[vol_pairwise_signif$vol_S2vsTd_pfdr<0.05,]

#Direction of results
vol_t_S2vsTd_neg <- is.negative(vol_t_S2vsTd_signif$vol_t_S2vsTd)


##Create table of region names and t-values for mapping in fslview
#Rename variables for merging
vol_t_S2vsTd_renamed <- rename(vol_t_S2vsTd_signif, ROI_NAME = volRegions)

#Remove the "mprage_jlf_vol_" prefix to match with the JLF labels
vol_t_S2vsTd_renamed[] <- lapply(vol_t_S2vsTd_renamed, function(x) gsub("mprage_jlf_vol_", "", x))

#Merge to add index numbers that correspond to the significant regions
vol_merged_S2vsTd <-merge(vol_t_S2vsTd_renamed,JLF.labels, by="ROI_NAME", all=FALSE)

#Remove the region names, leaving only the index numbers and t values
vol_merged_S2vsTd$ROI_NAME <- NULL

#Save as a .csv
write.table(vol_merged_S2vsTd, file="/data/jux/BBL/projects/pncHeterogeneity/subjectData/Tvalues/t_index_vol_S2vsTd.csv", row.names=FALSE, col.names=FALSE, quote = FALSE, sep=",")
