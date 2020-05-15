# TBMIBIquantifyPD1TcellsAcrossTopics.R
# Author: Erin McCaffrey
# Date created: 200225
# Overview: This script reads in the csv for the topic loadings for each cell with each cell assigned to a topic. 
# Next, it looks at the distribution of PD1+ T cells across topics in bulk.

library(ggplot2)
library(forcats)
library(dplyr)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Indicate cell type(s) to look at...##
cell_types<-c('lymphocyte')

##...Create dataframe with just cell type of interest...##
topic_data_cell <- topic_data[topic_data$cell_lin %in% cell_types, ]

##..Look at expression of marker on cell type across all topics..##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')
ggplot(data = topic_data_cell, aes(x = fct_reorder(as.factor(topic),PD.1,.fun=median,.desc=TRUE), y = PD.1, fill = as.factor(topic))) + 
  geom_violin(trim=FALSE, draw_quantiles = 0.5) + 
  labs(x = 'Topic') + 
  labs(y = 'PD1 Exprs') +
  scale_fill_manual(name = 'Topic', values = colors) 

##..Counts of all PD1 positive cells across topics..##

PD1_pos<-as.data.frame(table(topic_data_cell[topic_data_cell$PD.1>0.21, ]$topic))
colnames(PD1_pos)<-c('Topic','Count')

ggplot(data = PD1_pos, aes(x = as.factor(Topic), y = Count, fill = as.factor(Topic))) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'Count PD1+ Cells') +
  scale_fill_manual(name = 'Topic', values = colors) 

