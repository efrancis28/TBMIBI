# TBMIBIcalculateDiversityIndex.R
# Author: Erin McCaffrey 
# Date created: 190925
# Overview: This script reads in the immune cell frequency data, resizes the matrix to have
# each row be a sample and the columns be the absolute counts of all cell types (including fibroblasts,
# enothelium, and eptihelium). Next it calculates the Simpson and Shannon diversity index for each sample 
# and compares the trends in the sarcoid versus mycobacterial samples. 

library(factoextra)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)

# import data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

#cell frequencies
imm_freqs<-read.csv("immune_cell_freqs.csv")
lin_freqs<-read.csv("lineage_freqs.csv")

#reshape to a matrix where rows are samples and columns are counts of cell types
immfreq_data <- reshape2::dcast(imm_freqs, SampleID ~ cell_type, value.var = "count")
linfreq_data <- reshape2::dcast(lin_freqs, SampleID ~ lineage, value.var = "count")

#merge into one matrix with epithelial, endothelial, and fibroblast data
combo<-merge(immfreq_data,linfreq_data[,-5],by="SampleID")

#drop MAC

MAC_IDs<-c(57,58,19,87)
combo<-droplevels(combo[!combo$SampleID %in% MAC_IDs, ])

#add annotation of sample type
sarcoid<-c(67,68,69,70,71,72,73,74,75,76)

combo$Tissue <- "mycobacteria"
combo[combo$SampleID %in% sarcoid,]$Tissue <- "sarcoid"

#calculate Simpson diveristy index for each row
combo$simpsonIdx<-diversity(combo[,-c(1,22)], index = "simpson")
combo$shannonIdx<-diversity(combo[,-c(1,22,23)], index = "shannon")

#plot Simpson Index

color<-c("#9CD9D3","#00A79D")

simpson<- ggplot(combo, aes(x=Tissue, y=simpsonIdx, fill=Tissue)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  stat_compare_means(data=combo,mapping=aes(x=Tissue, y=simpsonIdx),label = "p.format", method= "wilcox.test",  comparisons = list(c("mycobacteria","sarcoid"))) +
  labs(x="Tissue Type") +
  labs(y="Simpson Index") 
simpson

#plot Simpson Index

shannon<- ggplot(combo, aes(x=Tissue, y=shannonIdx, fill=Tissue)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color)  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  stat_compare_means(data=combo,mapping=aes(x=Tissue, y=shannonIdx),label = "p.format", method= "t.test",  comparisons = list(c("mycobacteria","sarcoid"))) +
  labs(x="Tissue Type") +
  labs(y="Shannon Index") 
shannon

# save dataframe as csv

write.csv(combo, "sarcvTB_counts_diversity.csv", row.names = FALSE)


##supp
combo_freq<-combo
combo_freq[,2:21]<-combo_freq[,2:21]/rowSums(combo_freq[,-1])
combo_freq_TB<-droplevels(combo_freq[!combo_freq$SampleID %in% sarcoid,])
data_long <- gather(combo_freq_TB, factor_key=TRUE)

cell_stats<-data_long%>% group_by(key)%>%
  summarise(mean= mean(value), sd= sd(value), max = max(value),min = min(value))

