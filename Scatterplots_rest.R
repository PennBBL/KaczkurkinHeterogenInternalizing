###############################
#### REST ALFF SCATTERPLOT ####
###############################

#Load subject data
subjData <- readRDS("n840_DepAnxTd_rest_hydra_subjData.rds")

#Load the list of significant regions
signif_results <- read.csv("pfdr_omnibus_restAlff.csv", header=TRUE)

#Load libraries
library(mgcv)
library(mgcViz)

#The plot will jitter the points so they don't overlap as much, but we need to set a fixed seed for reproducibility
set.seed(123)

#Get significant variable names
restRegions <- as.vector(signif_results$restAlffRegions)

#Create a mean rest variable that includes all significant rest regions
subjData$averageRest_signifRegions <- rowMeans(subset(subjData, select = restRegions))

#Run model
restGam <- gam(averageRest_signifRegions ~ s(age) + sex + restRelMeanRMSMotion + Hydra_k2, method="REML", data=subjData)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
restGamViz <- getViz(restGam)

#Plot rest by group
plot(restGamViz, allTerms = TRUE, select = 4) + 
  labs(x = "", y = "Average ALFF (z)") + 
  l_fitPoints(mapping = aes(x=x,y=y,col=x)) + 
  l_points(mapping = aes(x=x,y=y,col=x), shape=19, size=2.5) + 
  l_ciBar(linetype=1, size=1) + ylim(-1, 1.75) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=33), 
        axis.text.x=element_text(size=30, colour="black"),
        axis.text.y=element_text(size=25, colour="black")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black"))

#Save scatterplot
ggsave(file="Scatterplot_rest.png")

