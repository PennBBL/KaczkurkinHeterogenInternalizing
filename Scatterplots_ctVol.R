################################
#### CT/VOLUME SCATTERPLOTS ####
################################

#Load data
subjData <- readRDS("n1141_DepAnxTd_t1_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(mgcViz)

#The plots will jitter the points so they don't overlap as much, but we need to set a fixed seed for reproducibility
set.seed(123)

##########
### CT ###
##########

##Get significant variable names (all ct regions were significant except L and R ENT)
ctRegions_all <- subjData[c(grep("mprage_jlf_ct",names(subjData)))]
ctRegions_signif <- ctRegions_all[,-grep("Ent",names(ctRegions_all))]
ctRegions <- names(ctRegions_signif)

#Create a mean ct variable that includes all significant ct regions
subjData$averageCT_signifRegions <- rowMeans(subset(subjData, select =ctRegions))

#Run model
ctGam <- gam(averageCT_signifRegions ~ s(age) + sex + averageManualRating + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
ctGamViz <- getViz(ctGam)

#Plot ct by group
plot(ctGamViz, allTerms = TRUE, select = 4) + 
  labs(x = "", y = "Average CT (mm)") +
  l_fitPoints(mapping = aes(x=x,y=y,col=x)) + 
  l_points(mapping = aes(x=x,y=y,col=x), shape=19, size=2.5) + 
  l_ciBar(linetype=1, size=1) + 
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=33), 
        axis.text.x=element_text(size=30, colour="black"),
        axis.text.y=element_text(size=25, colour="black")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black"))

#Save scatterplot
ggsave(file="Scatterplot_ct.png")

##############
### VOLUME ###
##############

#Divide total gray matter volume by 1000 to change the units from cubic millimeters (mm3) to cubic centimeters (cc3); 1 cc3 = 1,000 mm3
subjData$mprage_antsCT_vol_GrayMatter <- subjData$mprage_antsCT_vol_GrayMatter/1000

#Run model
volGam <- gam(mprage_antsCT_vol_GrayMatter ~ s(age) + sex + averageManualRating + Hydra_k2, method="REML", data=subjData)

##Convert the fitted object to the gamViz class to use the tools in mgcViz
volGamViz <- getViz(volGam)

##Plot total gray matter by group
plot(volGamViz, allTerms = TRUE, select = 4) + 
  labs(x = "", y = bquote('Total GM Volume'~(cm^3))) +
  l_fitPoints(mapping = aes(x=x,y=y,col=x)) + 
  l_points(mapping = aes(x=x,y=y,col=x), shape=19, size=2.5) + 
  l_ciBar(linetype=1, size=1) + 
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=33), 
        axis.text.x=element_text(size=30, colour="black"),
        axis.text.y=element_text(size=25, colour="black")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black"))

#Save scatterplot
ggsave(file="Scatterplot_vol.png")
