#TBMIBIexamineMETissueType.R
#Date created: 03/05/20
#Author: Erin McCaffrey
#This script reads in the topic data and assesses freqeucny of topics across three groups:
# 1. Lung resection 2. lung biopsy 3. extrapulmonary biopsy


library(ggpubr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(forcats)

##...Load in data..##

# Topic data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

# Define categories

# Plot individual examples and add origin data

resection<-c(64,65,21,84,42,88,28,89,85,13,35,36)
diagnostic_pulm<-c(6,7,14,21)
diagnostic_expulm<-c(33,34,26,27,40,61,47,48,54,55)

# Get topic frequencies

freq_topics <- as.data.frame(table(topic_data$SampleID, topic_data$topic))
colnames(freq_topics) <- c('SampleID', 'Topic', 'Count')
sample_totals<-as.data.frame(table(topic_data$SampleID))
freq_topics$Total <- rep(sample_totals$Freq, 8)
freq_topics$topic_freq <- freq_topics$Count / freq_topics$Total
freq_topics$Tissue<-'lung_resection'
freq_topics[freq_topics$SampleID %in% diagnostic_pulm,]$Tissue<-'lung_biopsy'
freq_topics[freq_topics$SampleID %in% diagnostic_expulm,]$Tissue<-'extrapulm_biopsy'

my_comparisons <- list(c("extrapulm_biopsy","lung_resection"),
                       c("extrapulm_biopsy","lung_biopsy"),
                       c("lung_biopsy","lung_resection"))


##..Produce plots with median line and dots..## 

data_summary <- freq_topics %>%
  group_by(Tissue,Topic) %>%
  summarize(combo_median = median(topic_freq),
            combo_se = sqrt(var(topic_freq)/length(topic_freq)))


ggplot(data = freq_topics, aes(x = Tissue, y = topic_freq)) + 
  geom_point(aes(x = Tissue, y = topic_freq), 
             position = position_jitter(width = 0.1, height = 0.0),size = 1) + 
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "red", width = 0.5, data = data_summary) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) + 
  stat_compare_means(label = "p.signif", method= "wilcox.test",comparisons=my_comparisons) +
  facet_wrap(~Topic, scale='free_y', ncol=4)
