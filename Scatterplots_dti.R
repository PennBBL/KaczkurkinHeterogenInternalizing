################################
#### DTI TRACTS SCATTERPLOT ####
################################

#Load subject data
subjData <- readRDS("n923_DepAnxTd_dti_hydra_subjData.rds")

#Load libraries
library(mgcv)
library(mgcViz)

#The plot will jitter the points so they don't overlap as much, but we need to set a fixed seed for reproducibility
set.seed(123)

#Get significant variable names
dtiRegions <- names(subjData)[grep("jhutract_fa_atr_|jhutract_fa_cst_r|jhutract_fa_cgh_|jhutract_fa_forceps_minor|jhutract_fa_ilf_|jhutract_fa_slf_r|jhutract_fa_uf_l",names(subjData))]

#Create a mean dti variable that includes all significant dti regions
subjData$averageDTI_signifRegions <- rowMeans(subset(subjData, select = dtiRegions))

#Run model
dtiGam <- gam(averageDTI_signifRegions ~ s(age) + sex + dti64Tsnr + Hydra_k2, method="REML", data=subjData)

#Run same model with Hydra_k2_reordered to make S2 the comparison group in order to see S1 vs S2 differences in the model.
dtiGam2 <- gam(averageDTI_signifRegions ~ s(age) + sex + dti64Tsnr + Hydra_k2_reordered, method="REML", data=subjData)

##Pass the results to the anova() function to fit ANOVAs (only anova() and not aov() can fit GAMs)
#This is your omnibus ANOVA test (tells you if the three groups are significantly different)
#df1=k-1, df2=n-k
dtiAnova <- anova(dtiGam)

##Pairwise comparisons
#Pull uncorrected p-values
dti_S1vsTd <- summary(dtiGam)$p.table[4,4]
dti_S2vsTd <- summary(dtiGam)$p.table[5,4]
dti_S1vsS2 <- summary(dtiGam2)$p.table[4,4]

#Combine the pairwise p values
dti_pairs <- cbind(dti_S1vsTd,dti_S2vsTd,dti_S1vsS2)

#Convert to data frame
dti_pairs <- as.data.frame(dti_pairs)

##FDR correction across the three pairwise group comparisons (S1vsTd, S2vsTd, and S1vsS2)
dti_fdr <- p.adjust(dti_pairs[1,],method="fdr")

#Print fdr-corrected p-values to three decimal places
dti_fdr_round <- round(dti_fdr,3)

#Convert the fitted object to the gamViz class to use the tools in mgcViz
dtiGamViz <- getViz(dtiGam)

#Plot dti by group
plot(dtiGamViz, allTerms = TRUE, select = 4) + 
  labs(x = "", y = "Average FA (z)") + 
  l_fitPoints(mapping = aes(x=x,y=y,col=x)) + 
  l_points(mapping = aes(x=x,y=y,col=x), shape=19, size=2.5) + 
  l_ciBar(linetype=1, size=1) +
  scale_y_continuous(limits = c(-2,2)) +
  scale_x_discrete(breaks=c("0","1","2"), labels=c("TD","S1","S2")) +
  theme(axis.title=element_text(size=33), 
        axis.text.x=element_text(size=30, colour="black"),
        axis.text.y=element_text(size=25, colour="black")) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), axis.line = element_line(colour="black"))

#Save scatterplot
ggsave(file="Scatterplot_dti.png")

