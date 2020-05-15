# TBMIBIdistributionofMyeloidCoreCellsAcrossMEs.R
# Author: Erin McCaffrey
# Date created: 200417
# Reads in the data for 18 FOVs with myeloid core annotation and the topic-annotated cell data. Next it looks at
# the frequency of cells in the myeloid core per topic.

library(ggplot2)
library(dplyr)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/summed_data")
region_data<-read.csv("celldata_region_annotated.csv")

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##..Subset just sample ID, cell label, and region..##

region_data_subset<-droplevels(region_data[,c("SampleID","cellLabelInImage","region")])

##..Subset the topic data to include on only samples with annotation..##

topic_data_subset<-droplevels(topic_data[topic_data$SampleID %in% unique(region_data$SampleID),])

##..Merge the region data..##

merged<-topic_data_subset %>% left_join(region_data_subset, by = c('SampleID','cellLabelInImage'))


##...Get frequency of cells in myeloid core per topic..##

#Across samples
freq_sample<-as.data.frame(table(merged$SampleID, merged$topic, merged$region))
total_sample<-as.data.frame(table(merged$SampleID, merged$topic))
freq_sample$topic_total<-rep(total_sample$Freq, n=4)
names(freq_sample)<-c("SampleID","topic","region","region_topic_count","topic_count")
freq_sample$freq_topic_region<-as.numeric(format((freq_sample$region_topic_count / freq_sample$topic_count)*100),digits=3)

##..Plot boxplot..##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

plot_data<-droplevels(freq_sample[freq_sample$region=='myeloid_core',])

ggplot(data = plot_data, aes(x = topic, y = freq_topic_region, fill = topic)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'ME') + 
  labs(y = 'Frequency of Cells in Myeloid Core') +
  scale_fill_manual(values = colors)

