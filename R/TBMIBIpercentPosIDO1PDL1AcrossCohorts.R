# TBMIBIpercentPosIDO1PDL1AcrossCohorts.R
# Author: Erin McCaffrey 
# Date created: 200317
# Overview: This script reads in the normalized data for IDO1 and PDL1. Next it plots the frequency of
# positive cells for IDO1 and PDL1 across the cohort broken down by sample type. 

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(ggpubr)

##..Import data, get just myeloid cells in TB and sarcoid..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")
MAC_IDs<-c(57,58,19,87)
sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
data<-droplevels(data[!data$SampleID %in% c(MAC_IDs,sarcoid), ])

##..Get the percent IDO and PDL1 pos cells in each point..#

IDO_thresh = 0.26
PDL1_thresh = 0.25

data_IDO<-droplevels(data[data$IDO>=IDO_thresh, ])
data_PDL1<-droplevels(data[data$PD.L1>=PDL1_thresh, ])


##..Across samples..##

freq_sample<-as.data.frame(table(data$SampleID))

#IDO
freq_IDO_sample<-as.data.frame(table(data_IDO$SampleID))
freq_sample$IDOpos<-freq_IDO_sample$Freq
names(freq_sample)<-c("SampleID","Total","Total_IDOpos")
freq_sample$percentIDO<-as.numeric(format((freq_sample$Total_IDOpos / freq_sample$Total)*100),digits=3)

#PDL1
freq_PDL1_sample<-as.data.frame(table(data_PDL1$SampleID))
freq_sample$Total_PDL1pos<-freq_PDL1_sample$Freq
freq_sample$percentPDL1<-as.numeric(format((freq_sample$Total_PDL1pos / freq_sample$Total)*100),digits=3)

##..Add sample type annotation..##

# Define categories

resection<-c(64,65,21,84,42,88,28,89,85,13,35,36)
diagnostic_pulm<-c(6,7,14,15)
diagnostic_expulm<-c(33,34,26,27,40,61,47,48,54,55)

freq_sample$Tissue<-'lung_resection'
freq_sample[freq_sample$SampleID %in% diagnostic_pulm,]$Tissue<-'lung_biopsy'
freq_sample[freq_sample$SampleID %in% diagnostic_expulm,]$Tissue<-'extrapulm_biopsy'

##..Plot..##

my_comparisons <- list(c("extrapulm_biopsy","lung_resection"),
                       c("extrapulm_biopsy","lung_biopsy"),
                       c("lung_biopsy","lung_resection"))

#IDO1

ggplot(data = freq_sample, aes(x = Tissue, y = percentIDO, fill = Tissue)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')


#PDL1

ggplot(data = freq_sample, aes(x = Tissue, y = percentPDL1, fill = Tissue)) + 
  geom_boxplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  labs(x = 'Tissue') + 
  labs(y = 'Frequency')

