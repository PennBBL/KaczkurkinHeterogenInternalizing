#################
### LOAD DATA ###
#################

##Demographic data (n=1629)
data.demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv", header=TRUE, na.strings="") 

##Clinical data
#Screening diagnoses (n=1601) (no missing values)
data.diag <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_psych_summary_vars_20131014.csv", header=TRUE)

#Psychosis clinical group (n=1601)
data.psychosis <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_diagnosis_dxpmr_20170509.csv", header=TRUE, na.strings="")

#Bifactors (n=1601)
data.bifactors <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_goassess_itemwise_bifactor_scores_20161219.csv", header=TRUE, na.strings="")

##Cognitive data
#Summary factor scores (n=1601)
data.cogFactors <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv", header=TRUE, na.strings="")

#WRAT scores (n=1601)
data.wrat <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_wrat_scores_20161215.csv", header=TRUE, na.strings="")

##Exclusion data
#Health exclusion (use the new healthExcludev2 variable) (n=1601; no missing values)
data.healthExclude <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20170421.csv", header=TRUE)

#T1 QA exclusion (n=1601)
data.t1QA <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv", header=TRUE, na.strings="NA")

#Resting state connectivity QA exclusion (n=1601)
data.restQA <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv", header=TRUE, na.strings="NA")

#DTI QA exclusion (n=1601)
data.dtiQA <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv", header=TRUE, na.strings="NA")

##Brain data
#JLF T1 ROIs (n=1601; no missing values)
data.ctJLF <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionCT_20170331.csv", header=TRUE)
data.volJLF <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_jlfAntsCTIntersectionVol_20170412.csv", header=TRUE)

#JLF total brain volume (TBV) (n=1601; no missing values)
data.tbv <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_ctVol20170412.csv", header=TRUE)

#JLF resting state connectivity ROIs (n=1601)
data.rest.alff <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_jlfALFFValues_20170714.csv", header=TRUE, na.strings="NA")

#JHU DTI ROIs (n=1601)
data.dti.FAtract <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_JHUTractFA_20170321.csv", header=TRUE, na.strings="NA")

#################
### DATA PREP ###
#################

#Transform the age variable from months to years
data.demo$age <- (data.demo$ageAtScan1)/12

#Define age squared (de-mean age)
data.demo$ageSq <- I(scale(data.demo$age, scale=FALSE, center=TRUE)^2)

#Recode male as 0 and female as 1 (0=male, 1=female)
data.demo$sex[which(data.demo$sex==1)] <- 0
data.demo$sex[which(data.demo$sex==2)] <- 1

#Make sex a factor
data.demo$sex <- as.factor(data.demo$sex)

#Define white vs nonwhite (white=1, non-white=0)
data.demo$white <- 0
data.demo$white[which(data.demo$race==1)] <- 1

#Make white a factor
data.demo$white <- as.factor(data.demo$white)

#z-score the resting state ALFF JLF regions and dti JHU tracts because they are on different scales
data.rest.alff[c(3:119)] <- lapply(data.rest.alff[c(3:119)], function(x) c(scale(x)))
data.dti.FAtract[c(3:20)] <- lapply(data.dti.FAtract[c(3:20)], function(x) c(scale(x)))

##################
### MERGE DATA ###
##################
dataMerge1 <-merge(data.demo,data.diag, by=c("bblid","scanid"), all=TRUE) 
dataMerge2 <-merge(dataMerge1,data.psychosis, by=c("bblid","scanid"), all=TRUE) 
dataMerge3 <-merge(dataMerge2,data.bifactors, by=c("bblid","scanid"), all=TRUE)
dataMerge4 <-merge(dataMerge3,data.cogFactors, by=c("bblid","scanid"), all=TRUE)
dataMerge5 <-merge(dataMerge4,data.wrat, by=c("bblid","scanid"), all=TRUE)
dataMerge6 <-merge(dataMerge5,data.healthExclude, by=c("bblid","scanid"), all=TRUE)
dataMerge7 <-merge(dataMerge6,data.t1QA, by=c("bblid","scanid"), all=TRUE)
dataMerge8 <-merge(dataMerge7,data.restQA, by=c("bblid","scanid"), all=TRUE)
dataMerge9 <- merge(dataMerge8,data.dtiQA, by=c("bblid","scanid"), all=TRUE)
dataMerge10 <- merge(dataMerge9,data.ctJLF, by=c("bblid","scanid"), all=TRUE)
dataMerge11 <- merge(dataMerge10,data.volJLF, by=c("bblid","scanid"), all=TRUE)
dataMerge12 <- merge(dataMerge11,data.tbv, by=c("bblid","scanid"), all=TRUE)
dataMerge13 <- merge(dataMerge12,data.rest.alff, by=c("bblid","scanid"), all=TRUE)
dataMerge14 <- merge(dataMerge13,data.dti.FAtract, by=c("bblid","scanid"), all=TRUE)

#Retain only the 1601 bblids (demographics has 1629)
data.n1601 <- dataMerge14[match(data.t1QA$bblid, dataMerge14$bblid, nomatch=0),] 

#Put bblids in ascending order
data.ordered <- data.n1601[order(data.n1601$bblid),]

#Count the number of subjects (should be 1601)
n <- nrow(data.ordered)

########################
### APPLY EXCLUSIONS ### 
########################
##Count the total number excluded for healthExcludev2=1 (1=Excludes those with medical rating 3/4, major incidental findings that distort anatomy, psychoactive medical medications)
#Included: n=1447; Excluded: n=154, but medical.exclude (n=81) + incidental.exclude (n=20) + medicalMed.exclude (n=64) = 165, so 11 people were excluded on the basis of two or more of these criteria
data.final <- data.ordered
data.final$ACROSS.INCLUDE.health <- 1
data.final$ACROSS.INCLUDE.health[data.final$healthExcludev2==1] <- 0
health.include<-sum(data.final$ACROSS.INCLUDE.health)
health.exclude<-1601-health.include

#Count the number excluded just medical rating 3/4 (GOAssess Medial History and CHOP EMR were used to define one summary rating for overall medical problems) (n=81)
data.final$ACROSS.INCLUDE.medical <- 1
data.final$ACROSS.INCLUDE.medical[data.final$medicalratingExclude==1] <- 0
medical.include<-sum(data.final$ACROSS.INCLUDE.medical)
medical.exclude<-1601-medical.include

#Count the number excluded for just major incidental findings that distort anatomy (n=20)
data.final$ACROSS.INCLUDE.incidental <- 1
data.final$ACROSS.INCLUDE.incidental[data.final$incidentalFindingExclude==1] <- 0
incidental.include<-sum(data.final$ACROSS.INCLUDE.incidental)
incidental.exclude<-1601-incidental.include

#Count the number excluded for just psychoactive medical medications (n=64)
data.final$ACROSS.INCLUDE.medicalMed <- 1
data.final$ACROSS.INCLUDE.medicalMed[data.final$psychoactiveMedMedicalv2==1] <- 0
medicalMed.include<-sum(data.final$ACROSS.INCLUDE.medicalMed)
medicalMed.exclude<-1601-medicalMed.include

#Subset the data to just those who pass healthExcludev2 (n=1447)
data.subset <-data.final[which(data.final$ACROSS.INCLUDE.health == 1), ]
n_health <- nrow(data.subset)

##Count the number excluded for failing to meet structural image quality assurance protocols
#Included: n=1396; Excluded: n=51
data.subset$ACROSS.INCLUDE.t1QA <- 1
data.subset$ACROSS.INCLUDE.t1QA[data.subset$t1Exclude==1] <- 0
t1QA.include<-sum(data.subset$ACROSS.INCLUDE.t1QA)
t1QA.exclude<-1447-t1QA.include

###Exclude those with ALL problems (health problems and problems with their t1 data) (included n=1396)
data.exclude <- data.subset[which(data.subset$healthExcludev2==0 & data.subset$t1Exclude == 0 ),]
n_health_t1 <- nrow(data.exclude)

##Count the number missing clinical data
#Included: n=1394; Excclinical.excludeluded: n=2
data.exclude$ACROSS.INCLUDE.clinical <- 1
data.exclude$ACROSS.INCLUDE.clinical[is.na(data.exclude$overall_psychopathology_4factorv2)] <- 0
clinical.include<-sum(data.exclude$ACROSS.INCLUDE.clinical)
clinical.exclude<-1396-clinical.include

#Exclude those missing clinical data
data.exclude <- data.exclude[!is.na(data.exclude$overall_psychopathology_4factorv2),]

#Check that number of subjects = 1394
n_clinical <- nrow(data.exclude)

##Count the number missing medu
#Included: n=1377; Excluded: n=17
data.exclude$ACROSS.INCLUDE.medu <- 1
data.exclude$ACROSS.INCLUDE.medu[is.na(data.exclude$medu1)] <- 0
medu.include<-sum(data.exclude$ACROSS.INCLUDE.medu)
medu.exclude<-1394-medu.include

#Exclude those missing medu
data.exclude <- data.exclude[!is.na(data.exclude$medu1),]

#Check that number of subjects = 1377
n_medu <- nrow(data.exclude)

##Count the number missing overall accuracy
#Included: n=1374; Excluded: n=3
data.exclude$ACROSS.INCLUDE.accuracy <- 1
data.exclude$ACROSS.INCLUDE.accuracy[is.na(data.exclude$Overall_Accuracy)] <- 0
accuracy.include<-sum(data.exclude$ACROSS.INCLUDE.accuracy)
accuracy.exclude<-1377-accuracy.include

#Exclude those missing accuracy data
data.exclude <- data.exclude[!is.na(data.exclude$Overall_Accuracy),]

#Check that number of subjects = 1374 
n_accuracy <- nrow(data.exclude)

##################################################
### DEFINE PSYCHOPATHOLOGY SCREENING DIAGNOSES ###
##################################################

##Make variables where 1 = diagnosis

subjData <- data.exclude

#ADHD
subjData$Add <- NA
subjData$Add[which(subjData$goassessSmryAdd==4)] <- 1

#Agoraphobia
subjData$Agr <- NA
subjData$Agr[which(subjData$goassessSmryAgr==4)] <- 1

#Anorexia
subjData$Ano <- NA
subjData$Ano[which(subjData$goassessSmryAno==4)] <- 1

#Bulimia
subjData$Bul <- NA
subjData$Bul[which(subjData$goassessSmryBul==4)] <- 1

#Conduct Disorder
subjData$Con <- NA
subjData$Con[which(subjData$goassessSmryCon==4)] <- 1

#Generalized Anxiety Disorder
subjData$Gad <- NA
subjData$Gad[which(subjData$goassessSmryGad==4)] <- 1

#Mania
subjData$Man <- NA
subjData$Man[which(subjData$goassessSmryMan==4)] <- 1

#Major Depressive Disorder
subjData$Mdd <- NA
subjData$Mdd[which(subjData$goassessSmryDep==4)] <- 1

#OCD
subjData$Ocd <- NA
subjData$Ocd[which(subjData$goassessSmryOcd==4)] <- 1

#Oppositional Defiant Disorder
subjData$Odd <- NA
subjData$Odd[which(subjData$goassessSmryOdd==4)] <- 1

#Panic Disorder
subjData$Pan <- NA
subjData$Pan[which(subjData$goassessSmryPan==4)] <- 1

#Psychosis
subjData$Ps <- NA
subjData$Ps[which(subjData$goassessDxpmr4=="4PS")] <- 1

#Posttraumatic Stress Disorder
subjData$Ptd <- NA
subjData$Ptd[which(subjData$goassessSmryPtd==4)] <- 1

#Separation Anxiety Disorder
subjData$Sep <- NA
subjData$Sep[which(subjData$goassessSmrySep==4)] <- 1

#Social Anxiety Disorder
subjData$Soc <- NA
subjData$Soc[which(subjData$goassessSmrySoc==4)] <- 1

#Specific Phobia
subjData$Sph <- NA
subjData$Sph[which(subjData$goassessSmryPhb==4)] <- 1

#Typically Developing
dxNames <- c("bblid","Add","Agr","Ano","Bul","Con","Gad","Man","Mdd","Ocd","Odd","Pan","Ps","Ptd","Sep","Soc","Sph")
dxDf <- data.matrix(subjData[,dxNames])
subjData$totDx <- rowSums(dxDf[,2:17], na.rm=TRUE) #This is how many people have how many diagnoses: sum(subjData$totDx==0): 428, sum(subjData$totDx==1): 321, sum(subjData$totDx>=2): 647.
subjData$Td <- 0
subjData$Td[which(subjData$totDx==0)] <- 1

#####################################
#### MAKE TD THE REFERENCE GROUP ####
#####################################

subjData$Add[which(subjData$Td==1)] <- 0
subjData$Agr[which(subjData$Td==1)] <- 0
subjData$Ano[which(subjData$Td==1)] <- 0
subjData$Bul[which(subjData$Td==1)] <- 0
subjData$Con[which(subjData$Td==1)] <- 0
subjData$Gad[which(subjData$Td==1)] <- 0
subjData$Man[which(subjData$Td==1)] <- 0
subjData$Mdd[which(subjData$Td==1)] <- 0
subjData$Ocd[which(subjData$Td==1)] <- 0
subjData$Odd[which(subjData$Td==1)] <- 0
subjData$Pan[which(subjData$Td==1)] <- 0
subjData$Ps[which(subjData$Td==1)] <- 0
subjData$Ptd[which(subjData$Td==1)] <- 0
subjData$Sep[which(subjData$Td==1)] <- 0
subjData$Soc[which(subjData$Td==1)] <- 0
subjData$Sph[which(subjData$Td==1)] <- 0

######################################################
### DEFINE DEPRESSION/ANXIETY VS TD GROUP VARIABLE ###
######################################################

##Create a variable (1=patient and 0=control) that will only include TD and those with Depression (Mdd) and/or Anx (Agr, Gad, Ocd, Pan, Ptd, Sep, Soc, Sph).
subjData$DepAnxTd <- NA
subjData$DepAnxTd[subjData$Agr==1|subjData$Gad==1|subjData$Mdd==1|subjData$Ocd==1|subjData$Pan==1|subjData$Ptd==1|subjData$Sep==1|subjData$Soc==1|subjData$Sph==1] <- 1
subjData$DepAnxTd[subjData$Td==1] <- 0

##Also create a DepAnxTd2 variable with 1=patient and -1=control for HYDRA.
subjData$DepAnxTd2 <- NA
subjData$DepAnxTd2[subjData$Agr==1|subjData$Gad==1|subjData$Mdd==1|subjData$Ocd==1|subjData$Pan==1|subjData$Ptd==1|subjData$Sep==1|subjData$Soc==1|subjData$Sph==1] <- 1
subjData$DepAnxTd2[subjData$Td==1] <- -1

#Make DepAnxTd and DepAnxTd2 factors
subjData$DepAnxTd <- as.factor(subjData$DepAnxTd)
subjData$DepAnxTd2 <- as.factor(subjData$DepAnxTd2)

#######################
### SAVE T1 DATASET ###
#######################

#Remove subjects with NAs, leaving only TD (0) and Anx/Dep (1)
subjData_t1 <- subjData[complete.cases(subjData$DepAnxTd),]

#Total sample size (n=1141)
n_t1 <- nrow(subjData_t1)

#See n's by group (DepAnx=715; TD=426)
ngrp_t1 <- table(subjData_t1$DepAnxTd)

#Save dataset
saveRDS(subjData_t1,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n1141_DepAnxTd_t1_subjData.rds")

###############################
#### SENSITIVITY EXCLUSION ####
###############################

##Count the number taking psychotropic psychiatric medications
#Included: n=1037; Excluded for meds: n=104
subjData_t1$ACROSS.INCLUDE.psychMeds <- 1
subjData_t1$ACROSS.INCLUDE.psychMeds[subjData_t1$psychoactiveMedPsychv2==1] <- 0
psychMeds.include<-sum(subjData_t1$ACROSS.INCLUDE.psychMeds)
psychMeds.exclude<-1141-psychMeds.include

#Subset the data to those not taking psychiatric meds
subjData_t1_noMeds <- subjData_t1[which(subjData_t1$ACROSS.INCLUDE.psychMeds==1),]

###################################
### SAVE T1 SENSITIVITY DATASET ###
###################################

#Total sample size (n=1037)
n_t1_noMeds <- nrow(subjData_t1_noMeds)

#See n's by group (DepAnx=628; TD=409)
ngrp_t1_noMeds <- table(subjData_t1_noMeds$DepAnxTd)

#Save dataset
saveRDS(subjData_t1_noMeds,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n1037_DepAnxTd_t1_noMeds_subjData.rds")

########################
#### REST EXCLUSION ####
########################

#Resting state connectivity has the additional exclusion criteria: restExclude.
#Count the number excluded for failing to meet resting state image quality assurance protocols
#Included: n=840; Excluded: n=301
subjData_t1$ACROSS.INCLUDE.restQA <- 1
subjData_t1$ACROSS.INCLUDE.restQA[subjData_t1$restExclude==1] <- 0
restQA.include<-sum(subjData_t1$ACROSS.INCLUDE.restQA)
restQA.exclude<-1141-restQA.include

#Subset the data to just those who pass restQA
subjData_rest <- subjData_t1[which(subjData_t1$ACROSS.INCLUDE.restQA==1),]

#########################
### SAVE REST DATASET ###
#########################

#Total sample size (n=840)
n_rest <- nrow(subjData_rest)

#See n's by group (DepAnx=530; TD=310)
ngrp_rest <- table(subjData_rest$DepAnxTd)

#Save dataset
saveRDS(subjData_rest,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n840_DepAnxTd_rest_subjData.rds")

#####################################
### SAVE REST SENSITIVITY DATASET ###
#####################################

#Subset the data to those not taking psychiatric meds (excluded: n=77)
subjData_rest_noMeds <- subjData_rest[which(subjData_rest$ACROSS.INCLUDE.psychMeds==1),]

#Total sample size (n=763)
n_rest_noMeds <- nrow(subjData_rest_noMeds)

#See n's by group (DepAnx=466; TD=297)
ngrp_rest_noMeds <- table(subjData_rest_noMeds$DepAnxTd)

#Save dataset
saveRDS(subjData_rest_noMeds,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n763_DepAnxTd_rest_noMeds_subjData.rds")

#######################
#### DTI EXCLUSION ####
#######################

#DTI has the additional exclusion criteria: dti64Exclude.
#Count the number excluded for failing to meet DTI image quality assurance protocols
#Included: n=923; Excluded: n=218
subjData_t1$ACROSS.INCLUDE.dtiQA <- 1
subjData_t1$ACROSS.INCLUDE.dtiQA[subjData_t1$dti64Exclude==1] <- 0
dtiQA.include<-sum(subjData_t1$ACROSS.INCLUDE.dtiQA)
dtiQA.exclude<-1141-dtiQA.include

#Subset the data to just those who pass dtiQA
subjData_dti <- subjData_t1[which(subjData_t1$ACROSS.INCLUDE.dtiQA==1),]

########################
### SAVE DTI DATASET ###
########################

#Total sample size (n=923)
n_dti <- nrow(subjData_dti)

#See n's by group (DepAnx=579; TD=344)
ngrp_dti <- table(subjData_dti$DepAnxTd)

#Save dataset
saveRDS(subjData_dti,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n923_DepAnxTd_dti_subjData.rds")

####################################
### SAVE DTI SENSITIVITY DATASET ###
####################################

#Subset the data to those not taking psychiatric meds (excluded: n=87)
subjData_dti_noMeds <- subjData_dti[which(subjData_dti$ACROSS.INCLUDE.psychMeds==1),]

#Total sample size (n=836)
n_dti_noMeds <- nrow(subjData_dti_noMeds)

#See n's by group (DepAnx=507; TD=329)
ngrp_dti_noMeds <- table(subjData_dti_noMeds$DepAnxTd)

#Save dataset
saveRDS(subjData_dti_noMeds,"/data/jux/BBL/projects/pncHeterogeneity/subjectData/n836_DepAnxTd_dti_noMeds_subjData.rds")
