############################################################
#### GROUP DIFFERENCES IN COGNITION AND PSYCHOPATHOLOGY ####
############################################################

#Load data
subjData <- readRDS("n1141_DepAnxTd_t1_hydra_subjData.rds")

#Check that n=1141
n <- nrow(subjData)

#Load libraries
library(mgcv)
library(mgcViz)

#The plot will jitter the points so they don't overlap as much, but we need to set a fixed seed for reproducibility
set.seed(123)

#################################################
### EXECUTIVE FUNCTIONING COGNITIVE REASONING ###
#################################################

#Run model
execLm <- gam(F1_Exec_Comp_Res_Accuracy ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
execLmViz <- getViz(execLm)

#Plot executive functioning by group
plot(execLmViz, allTerms = TRUE, select = 3) +
  ylab("Ex Function Estimates") + xlab("") + ylim(-.5, .3) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_ExecFuncCogReasoning.png")

########################
### SOCIAL COGNITION ###
########################

#Run model
socLm <- gam(F2_Social_Cog_Accuracy ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
socLmViz <- getViz(socLm)

#Plot social cognition by group
plot(socLmViz, allTerms = TRUE, select = 3) +
  ylab("Social Cog Estimates") + xlab("") + ylim(-.4, .4) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_SocialCognition.png")

#######################
### EPISODIC MEMORY ###
#######################

#Run model
epiLm <- gam(F3_Memory_Accuracy ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
epiLmViz <- getViz(epiLm)

#Plot episodic memory by group
plot(epiLmViz, allTerms = TRUE, select = 3) +
  ylab("Memory Estimates") + xlab("") + ylim(-.4, .4) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_EpisodicMemory.png")

###############################
### OVERALL PSYCHOPATHOLOGY ###
###############################

#Run model
psychLm <- gam(overall_psychopathology_4factorv2 ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
psychLmViz <- getViz(psychLm)

#Plot overall psychopathology by group
plot(psychLmViz, allTerms = TRUE, select = 3) +
  ylab("Overall Psy Estimates") + xlab("") + 
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_overallPsychopathology.png")

##################
### BEHAVIORAL ###
##################

#Run model
extLm <- gam(externalizing_4factorv2 ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
extLmViz <- getViz(extLm)

#Plot externalizing by group
plot(extLmViz, allTerms = TRUE, select = 3) +
  ylab("Behavioral Estimates") + xlab("") + ylim(-.25, .8) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_externalizing.png")

############
### FEAR ###
############

#Run model
fearLm <- gam(phobias_4factorv2 ~ age + sex + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
fearLmViz <- getViz(fearLm)

#Plot fear by group
plot(fearLmViz, allTerms = TRUE, select = 3) +
  ylab("Fear Estimates") + xlab("") + ylim(-.25, .8) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x), shape=19, size=7) +
  l_ciBar(linetype=1, size=2, col=c("#f35e5a","#17b12b","#5086ff")) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=40),
        axis.text.x=element_text(size=35, face="bold", colour=c("#f35e5a","#17b12b","#5086ff")),
        axis.text.y=element_text(size=30, face="bold", colour="black"),
        axis.title.y = element_text(face="bold")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(size=2, colour="black"))

#Save plot
ggsave(file="GrpDiff_fear.png")
