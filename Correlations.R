######################
#### Correlations ####
######################

#Load data
subjData <- readRDS("/data/jux/BBL/projects/pncHeterogeneity/subjectData/MergedWithHYDRA/n1141_DepAnxTd_t1_hydra_ICV_subjData.rds")

#Load libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(ggcorrplot)

########################
### DEFINE VARIABLES ###
########################

#Get CT variables
ctRegions <- names(subjData)[grep("mprage_jlf_ct",names(subjData))]

#Create averageCT
subjData$averageCT <- rowMeans(subset(subjData, select = ctRegions))

#Get ALFF variables
restAlff_all <- subjData[c(grep("rest_jlf_alff",names(subjData)))]
restAlff_short <- restAlff_all[,-grep("Cerebell",names(restAlff_all))]
restAlffRegions <- names(restAlff_short)

#Create	averageAlff
subjData$averageAlff <- rowMeans(subset(subjData, select = restAlffRegions))

#Get FA variables
dtiTrRegions <- names(subjData)[grep("dti_dtitk_jhutract",names(subjData))]

#Create averageFA
subjData$averageFA <- rowMeans(subset(subjData, select = dtiTrRegions))

####################
### DESCRIPTIVES ###
####################

#Means and SDs
meanSdICV <- ddply(subjData,~Hydra_k2,summarise,mean=mean(ICV),sd=sd(ICV))
meanSdTBV <- ddply(subjData,~Hydra_k2,summarise,mean=mean(mprage_antsCT_vol_TBV),sd=sd(mprage_antsCT_vol_TBV))
meanSdTotGray <- ddply(subjData,~Hydra_k2,summarise,mean=mean(mprage_antsCT_vol_GrayMatter),sd=sd(mprage_antsCT_vol_GrayMatter))
meanSdAvgCT <- ddply(subjData,~Hydra_k2,summarise,mean=mean(averageCT),sd=sd(averageCT))
meanSdAvgAlff <- ddply(subjData,~Hydra_k2,summarise,mean=mean(averageAlff,na.rm=TRUE),sd=sd(averageAlff,na.rm=TRUE))
meanSdAvgFA <- ddply(subjData,~Hydra_k2,summarise,mean=mean(averageFA,na.rm=TRUE),sd=sd(averageFA,na.rm=TRUE))

####################
### CORRELATIONS ###
####################

#Correlation between ICV and TBV
corr_ICV_TBV <- cor(subjData$ICV, subjData$mprage_antsCT_vol_TBV)

#Correlation between ICV and total gray matter volume
corr_ICV_totGray <- cor(subjData$ICV, subjData$mprage_antsCT_vol_GrayMatter)

############################
### CORRELATION MATRICES ###
############################

##Subset only by variables needed for correlation tables
VarNames <- c("Hydra_k2","ICV","mprage_antsCT_vol_TBV","mprage_antsCT_vol_GrayMatter","averageCT","averageAlff","averageFA","Overall_Accuracy","F1_Exec_Comp_Res_Accuracy","F2_Social_Cog_Accuracy","F3_Memory_Accuracy")
Vars <- subjData[VarNames]

#Rename variables for plotting
names(Vars)[names(Vars) == 'mprage_antsCT_vol_TBV'] <- 'TBV'
names(Vars)[names(Vars) == 'mprage_antsCT_vol_GrayMatter'] <- 'TotalGrayMatter'
names(Vars)[names(Vars) == 'averageCT'] <- 'CT'
names(Vars)[names(Vars) == 'averageAlff'] <- 'ALFF'
names(Vars)[names(Vars) == 'averageFA'] <- 'FA'
names(Vars)[names(Vars) == 'Overall_Accuracy'] <- 'CogAccuracy'
names(Vars)[names(Vars) == 'F1_Exec_Comp_Res_Accuracy'] <- 'ExFunction'
names(Vars)[names(Vars) == 'F2_Social_Cog_Accuracy'] <- 'SocialCog'
names(Vars)[names(Vars) == 'F3_Memory_Accuracy'] <- 'Memory'

#Split file by group to get separate correlation tables for S1, S2, and Td
#S1
VarsS1 <- Vars[which(Vars$Hydra_k2==1),]

#Check that n=403
nVarsS1 <- nrow(VarsS1)

#Remove Hydra_k2 variable
VarsS1$Hydra_k2 <- NULL

#S2
VarsS2 <- Vars[which(Vars$Hydra_k2==2),]

#Check that n=312
nVarsS2 <- nrow(VarsS2)

#Remove Hydra_k2 variable
VarsS2$Hydra_k2 <- NULL

#TD
VarsTd <- Vars[which(Vars$Hydra_k2==0),]

#Check that n=426
nVarsTd <- nrow(VarsTd)

#Remove Hydra_k2 variable
VarsTd$Hydra_k2 <- NULL

##S1 corrplot
#Compute a correlation matrix
corr_S1 <- round(cor(VarsS1, use="complete.obs"), 2)

#Compute a matrix of correlation p-values
p_S1 <- cor_pmat(VarsS1)

#Plot heat map with correlation values
ggcorrplot(corr_S1, p.mat = p_S1, type = "lower", insig = "blank", lab = TRUE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/Corrs_S1.png")

##S2 corrplot
#Compute a correlation matrix
corr_S2 <- round(cor(VarsS2, use="complete.obs"), 2)

#Compute a matrix of correlation p-values
p_S2 <- cor_pmat(VarsS2)

#Plot heat map with correlation values
ggcorrplot(corr_S2, p.mat = p_S2, type = "lower", insig = "blank", lab = TRUE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/Corrs_S2.png")

##TD corrplot
#Compute a correlation matrix
corr_Td <- round(cor(VarsTd, use="complete.obs"), 2)

#Compute a matrix of correlation p-values
p_Td <- cor_pmat(VarsTd)

#Plot heat map with correlation values
ggcorrplot(corr_Td, p.mat = p_Td, type = "lower", insig = "blank", lab = TRUE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 
ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/Corrs_Td.png")

##Subset to only ICV, TBV, total gray matter, CT, ALFF, and FA
VarNames2 <- c("ICV","TBV","TotalGrayMatter","CT","ALFF","FA")
Vars2 <- Vars[VarNames2]

##Corrplot (across all groups)
#Compute a correlation matrix
corr <- round(cor(Vars2, use="complete.obs"), 2)

#Compute a matrix of correlation p-values
p <- cor_pmat(Vars2)

#Plot heat map with correlation values
ggcorrplot(corr, p.mat = p, type = "lower", insig = "blank", lab = TRUE) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/Corrs_AllGroups.png")


####################
### SCATTERPLOTS ###
####################

#Scatterplot of ICV by TBV
ggplot(subjData, aes(x=mprage_antsCT_vol_TBV, y=ICV)) + geom_point() + 
  theme_classic(base_size = 20) + labs(x = "TBV") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))

ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/ICVbyTBV.png")

#Scatterplot of ICV by total gray matter volume
ggplot(subjData, aes(x=mprage_antsCT_vol_GrayMatter, y=ICV)) + geom_point() +
  theme_classic(base_size = 20) + labs(x = "Total Gray Matter Volume") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))	+
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE))

ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/ICVbyTotGray.png")
