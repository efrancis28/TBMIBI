# TBMIBIquantifyCellAbundanceinMEPerPatient.R
# Author: Erin McCaffrey
# Date created: 191205
# Overview: This script reads in the csv for the topic loadings for each cell. It assigns each cell to a topic 
# (max loading) and then plots the abundance of each cell type in each topic across each patient (ie. one subplot per patient and topic)

library(ggplot2)
library(dplyr)

##...Load in data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Quantify cell type abundances in each topic...##

topic_freqs<-as.data.frame(table(topic_data$cell_type, topic_data$topic, topic_data$SampleID))
colnames(topic_freqs)<-c('cell_type','topic','SampleID','count')

##...Plot broken down by cell type..##

imm_order <- unique(topic_freqs$cell_type)
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
colorkey<-read.csv('colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% imm_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = imm_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

point_order_gran<-c(64, 65, 21, 84, 42, 88, 28, 89, 85, 13, 35, 36, 14, 15, 57, 58, 19, 87, 6, 7, 33, 34, 26, 27, 40, 61, 47, 48, 54, 55)
topic_freqs$SampleID <- factor(topic_freqs$SampleID, levels=point_order_gran)

ggplot(data = topic_freqs, aes(x = SampleID, y = count, fill = cell_type)) + 
  geom_col() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x = 'SampleID') + 
  labs(y = 'Cell Count') +
  scale_fill_manual(name = 'Topic', values = color) +
  facet_wrap(.~topic, scale='free_y', ncol = 2)

