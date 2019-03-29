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

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
ctGam2 <- gam(averageCT_signifRegions ~ s(age) + sex + averageManualRating + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
ctAnova <- anova(ctGam)

##Pairwise comparisons
#Pull uncorrected p-values
ct_S1vsTd <- summary(ctGam)$p.table[4,4]
ct_S2vsTd <- summary(ctGam)$p.table[5,4]
ct_S1vsS2 <- summary(ctGam2)$p.table[4,4]

#Combine the pairwise p values
ct_pairs <- cbind(ct_S1vsTd,ct_S2vsTd,ct_S1vsS2)

#Convert to data frame
ct_pairs <- as.data.frame(ct_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
ct_fdr <- p.adjust(ct_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
ct_fdr_round <- round(ct_fdr,3)

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

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
volGam2 <- gam(mprage_antsCT_vol_GrayMatter ~ s(age) + sex + averageManualRating + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
volAnova <- anova(volGam)

##Pairwise comparisons
#Pull uncorrected p-values
vol_S1vsTd <- summary(volGam)$p.table[4,4]
vol_S2vsTd <- summary(volGam)$p.table[5,4]
vol_S1vsS2 <- summary(volGam2)$p.table[4,4]

#Combine the pairwise p values
vol_pairs <- cbind(vol_S1vsTd,vol_S2vsTd,vol_S1vsS2)

#Convert to data frame
vol_pairs <- as.data.frame(vol_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
vol_fdr <- p.adjust(vol_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
vol_fdr_round <- round(vol_fdr,3)

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
