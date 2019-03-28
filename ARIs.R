###################
#### PLOT ARIs ####
###################

#Load data (NOTE: This is a modified version of DepAnxTd_ctVolJLF_ARI.csv with headers (subtype and ARI), subtype numbers added, and subtype 1=0 removed because HYDRA does not compute an ARI for 1 subtype)
subjData <- read.csv("/data/jux/BBL/projects/pncHeterogeneity/subjectData/ctVolJLF_ARIs.csv", header=TRUE)

#Load libraries
library(ggplot2)

ggplot(data=subjData, aes(x=subtype, y=ARI, group=1)) +
  geom_line(color="black", size=1.5) +
  theme_classic(base_size = 24) +
  theme(legend.position="none") +
  scale_y_continuous(breaks=c(0,.10,.20,.30,.40,.50,.60,.70), 
  labels = scales::number_format(accuracy = 0.01), expand = c(0, 0)) +
  expand_limits(y=.75) +
  scale_x_continuous(breaks=c(2,3,4,5,6,7,8,9,10)) +
  labs(x = "Number of Subtypes") +
  geom_segment(aes(x = 2, y = 0, xend = 2, yend = .66), linetype="dashed", color = "red", size=1.5) +
  annotate(geom="text", x=2.2, y=.69, label="0.66", color="red", size=7.5) 

ggsave(file="/data/jux/BBL/projects/pncHeterogeneity/tablesFigures/ARIs.png")
