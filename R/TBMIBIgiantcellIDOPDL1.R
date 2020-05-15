# TBMIBIgiantcellIDOPDL1.R
# Author: Erin McCaffrey 
# Date created: 190925
# Overview: Reads in the normalized intensity and annotated dataframe. Determines the frequency of 
# IDO and PDL1 pos cells in giant cells. 

library(ggplot2)

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

data_norm<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")

##..Get percent IDO and PDL1 positive.##

IDO_thresh = 0
PDL1_thresh = 0

##..Get total number of IDO1 and PDL1+ giant cells..##

GC<-as.data.frame(as.numeric(sum(data_norm$cell_type=="giant_cell")))
colnames(GC)<-"GC_total"
GC$IDOpos<-as.numeric(sum(data_norm[data_norm$cell_type=="giant_cell",]$IDO>IDO_thresh))
GC$PDL1pos<-as.numeric(sum(data_norm[data_norm$cell_type=="giant_cell",]$PD.L1>PDL1_thresh))

##..Get percent IDO1 and PDL1 positive..##

GC$IDOpos_percent<-GC$IDOpos/GC$GC_total
GC$PDL1pos_percent<-GC$PDL1pos/GC$GC_total

markers<-c("IDO","PDL1")
percent<-c(GC$IDOpos_percent,GC$PDL1pos_percent)
GC_markers<-data.frame(markers,percent)

##..Plot..##

bulk<-ggplot(data=GC_markers, aes(x=markers,y=percent)) +
  geom_bar(stat="identity", fill="grey", lwd=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="Marker") +
  labs(y="% Giant Cells Positive") +
  guides(fill=FALSE) 
bulk
dev.off()