# TBMIBIquantifyLymphMEEnrichments.R
# Author: Erin McCaffrey
# Date created: 200109
# Overview: This script reads in the csv for the topic loadings for each cell with each cell assigned to a topic. 
# Next, for a given marker and cell population it produces a violin plot of expression for that marker by cells across
# all topics to assess cell-intrinsic expression profiles versus niche-specific expression profiles.

library(ggpubr)
library(ggplot2)
library(dplyr)
library(forcats)
library(splitstackshape)
library(reshape2)


##..Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##..Filter non-lymphocyte cell types and also an immune cell only version..##

topic_data_imm <-  droplevels(topic_data[topic_data$lineage == 'immune', ])
cell_types<-c('lymphocyte')
topic_data <- droplevels(topic_data[topic_data$cell_lin %in% cell_types, ])


##..Get count of each cell-type as a function of topic..##

topic_freqs_lymph<-as.data.frame(table(topic_data$SampleID,topic_data$cell_type,topic_data$topic))
colnames(topic_freqs_lymph) <- c('SampleID','cell_type','topic','count')

##..Get frequency of each cell type out of total for that cell type per ROI..##

cell_totals<-as.data.frame(table(topic_data$SampleID,topic_data$cell_type))
totals<-rep(cell_totals$Freq, 8)
topic_freqs_lymph$cell_type_total<-totals
topic_freqs_lymph$cell_type_freq <- as.numeric(format((topic_freqs_lymph$count / topic_freqs_lymph$cell_type_total)),digits=3)

##..Get frequency the frequency of immune out of total immune for each topic per ROI..##

imm_freqs <- as.data.frame(table(topic_data_imm$SampleID,topic_data_imm$topic))
colnames(imm_freqs) <- c('SampleID','topic','imm_count')
imm_totals<-as.data.frame(table(topic_data_imm$SampleID))
imm_total_counts<-rep(imm_totals$Freq,8)
imm_freqs$total_imm<-imm_total_counts
imm_freqs$freq_of_imm <- as.numeric(format((imm_freqs$imm_count / imm_freqs$total_imm)),digits=3)

##..Get frequency the frequency of lymphs out of total lymphs for each topic per ROI..##

lymph_freqs <- as.data.frame(table(topic_data$SampleID,topic_data$topic))
colnames(lymph_freqs) <- c('SampleID','topic','lymph_count')
lymph_totals<-as.data.frame(table(topic_data$SampleID))
lymph_total_counts<-rep(lymph_totals$Freq,8)
lymph_freqs$total_lymph<-lymph_total_counts
lymph_freqs$freq_of_lymph <- as.numeric(format((lymph_freqs$lymph_count / lymph_freqs$total_lymph)),digits=3)

##..Get frequency of T cells out of total T cells for each topic..##

topic_data_T <-droplevels(topic_data[topic_data$cell_type != 'B_cell',])
T_freqs <- as.data.frame(table(topic_data_T$SampleID, topic_data_T$topic))
colnames(T_freqs) <- c('SampleID','topic','T_count')
T_totals<-as.data.frame(table(topic_data_T$SampleID))
T_total_counts<-rep(T_totals$Freq,8)
T_freqs$total_T<-T_total_counts
T_freqs$freq_of_T <- as.numeric(format((T_freqs$T_count/ T_freqs$total_T)),digits=3)

##..Get frequency of CD4 T cells out of total CD4 T cells for each topic..##

topic_data_CD4 <-droplevels(topic_data[topic_data$cell_type == 'CD4_T',])
CD4T_freqs <- as.data.frame(table(topic_data_CD4$SampleID, topic_data_CD4$topic))
colnames(CD4T_freqs) <- c('SampleID','topic','CD4T_count')
CD4T_totals<-as.data.frame(table(topic_data_CD4$SampleID))
CD4T_total_counts<-rep(CD4T_totals$Freq,8)
CD4T_freqs$total_CD4<-CD4T_total_counts
CD4T_freqs$freq_of_CD4 <- as.numeric(format((CD4T_freqs$CD4T_count/ CD4T_freqs$total_CD4)),digits=3)

##..Reorder the summary sheet to append the various baseline data..##

topic_freqs_lymph<-topic_freqs_lymph[order(topic_freqs_lymph$SampleID),]

##..Get log 2 FC of cell type frequency of total cell type to frequency of total immune..##

imm_freqs  <- expandRows(imm_freqs, count = 4, count.is.col = F, drop = F)
imm_freqs <-imm_freqs [order(imm_freqs $SampleID),]
topic_freqs_lymph$total_imm_freq <- imm_freqs$freq_of_imm
topic_freqs_lymph$total_imm_FC <- log2(topic_freqs_lymph$cell_type_freq / topic_freqs_lymph$total_imm_freq)

##..Get log 2 FC of cell type frequency of total cell type to frequency of total lymphs..##

lymph_freqs <- expandRows(lymph_freqs, count = 4, count.is.col = F, drop = F)
lymph_freqs <-lymph_freqs[order(lymph_freqs$SampleID),]
topic_freqs_lymph$total_lymph_freq <- lymph_freqs$freq_of_lymph
topic_freqs_lymph$total_lymph_FC <- log2(topic_freqs_lymph$cell_type_freq / topic_freqs_lymph$total_lymph_freq)

##..Get log 2 FC of cell type frequency of total cell type to frequency of total T..##

T_freqs <- expandRows(T_freqs, count = 4, count.is.col = F, drop = F)
T_freqs <-T_freqs[order(T_freqs$SampleID),]
topic_freqs_lymph$total_T_freq <- T_freqs$freq_of_T
topic_freqs_lymph$total_T_FC <- log2(topic_freqs_lymph$cell_type_freq / topic_freqs_lymph$total_T_freq)

##..Get log 2 FC of cell type frequency of total cell type to frequency of total T..##

CD4T_freqs <- expandRows(CD4T_freqs, count = 4, count.is.col = F, drop = F)
CD4T_freqs <-CD4T_freqs[order(CD4T_freqs$SampleID),]
topic_freqs_lymph$total_CD4T_freq <- CD4T_freqs$freq_of_CD4
topic_freqs_lymph$total_CD4T_FC <- log2(topic_freqs_lymph$cell_type_freq / topic_freqs_lymph$total_CD4T_freq)

##..Plot the FC relative to all immune, lymphs, T cells across topics for each cell type..##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

ggplot(data = topic_freqs_lymph, aes(x = cell_type, y = total_imm_FC, fill = topic)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors) +
  labs(x = 'Topic') + 
  labs(y = 'log2(Frequency of Subset/ Frequency of Total Immune)') +
  facet_wrap(~topic, scale='free_y')

plot_colors<-c('#CA6627', '#7470AE', '#D53E88', '#74A439')
plot_topics<-c(1,2,3,4)
plot_data<-droplevels(topic_freqs_lymph[topic_freqs_lymph$topic %in% plot_topics,])

# Plot CD4 and CD8 T cells relative to total immune

CD8_plot<-plot_data %>% filter((!is.infinite(total_imm_FC)))
CD8_plot<-CD8_plot[CD8_plot$cell_type=='CD8_T',]

ggplot(data = CD8_plot, aes(x = topic, y = total_imm_FC, fill = topic)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = plot_colors) +
  geom_point(size = 3, position = position_jitterdodge()) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'log2(Frequency of Subset/ Frequency of Total Immune)')

CD4_plot<-plot_data %>% filter((!is.infinite(total_imm_FC)))
CD4_plot<-CD4_plot[CD4_plot$cell_type=='CD4_T',]

ggplot(data = CD4_plot, aes(x = topic, y = total_imm_FC, fill = topic)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = plot_colors) +
  geom_point(size = 3, position = position_jitterdodge()) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'log2(Frequency of Subset/ Frequency of Total Immune)')


# Plot Tregs relative to CD4 T cells

treg_plot<-plot_data %>% filter((!is.infinite(total_CD4T_FC)))

ggplot(data = treg_plot, aes(x = topic, y = total_CD4T_FC, fill = topic)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = plot_colors) +
  geom_point(size = 3, position = position_jitterdodge()) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  stat_compare_means(label = "p.format", method= "wilcox.test",comparisons=list(c(1,3))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'log2(Frequency of Subset/ Frequency of Total CD4')

  

  
