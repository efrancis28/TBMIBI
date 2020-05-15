# TBMIBIquantifyPDL1PosFreqAcrossMEsAndCells.R
# Author: Erin McCaffrey
# Date created: 200106

library(ggpubr)
library(ggplot2)
library(colorspace)
library(forcats)
library(dplyr)
library(tibble)
library(reshape2)
library(gplots)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Filter non-myeloid cell types...##

cell_types<-c('myeloid')
topic_data <- droplevels(topic_data[topic_data$cell_lin %in% cell_types, ])


##...Get percent of cells positive for PDL1 and ratio of pos:neg...##

marker_thresh <- 0.25
data_marker<-droplevels(topic_data[topic_data$PD.L1>marker_thresh, ])

#Across samples
freq_sample<-as.data.frame(table(topic_data$SampleID, topic_data$topic, topic_data$cell_type))
freq_marker_sample<-as.data.frame(table(data_marker$SampleID, data_marker$topic, data_marker$cell_type))
freq_sample$PD.L1pos <- freq_marker_sample$Freq
names(freq_sample)<-c("SampleID","topic","cell_type","Total","Total_PDL1pos")
freq_sample$Total_PDL1neg <- freq_sample$Total - freq_sample$Total_PDL1pos
freq_sample$percentPDL1pos<-as.numeric(format((freq_sample$Total_PDL1pos / freq_sample$Total)*100),digits=3)
freq_sample$FC <- log2(freq_sample$Total_PDL1pos/freq_sample$Total_PDL1neg)

##...Visualize plot of FC across topics (all cell types, boxplot...##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

ggplot(data = freq_sample, aes(x = topic, y = FC, fill = topic)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'log2(PDL1+ / PDL1-)') +
  scale_fill_manual(values = colors)

##...Produce table of frequency of PDL1+ for all myeloid cells and then each topic...##

# across topics and points
freq_topic<-as.data.frame(table(topic_data$topic, topic_data$SampleID))
freq_marker_topic<-as.data.frame(table(data_marker$topic, data_marker$SampleID))
freq_topic$PDL1pos <- freq_marker_topic$Freq
names(freq_topic)<-c("topic","SampleID","Total","Total_PDL1pos")

# total myeloid
freq_myeloid<-as.data.frame(table(topic_data$SampleID))
freq_myeloid<-add_column(freq_myeloid, topic = 'total_myeloid', .before = "Var1")
freq_myeloid_marker<-as.data.frame(table(data_marker$SampleID))
freq_myeloid$PDL1pos <- freq_myeloid_marker$Freq
names(freq_myeloid)<-c("topic","SampleID","Total","Total_PDL1pos")

# merge and get frequency
freq_topic<-rbind(freq_myeloid, freq_topic)
freq_topic$percentPDL1pos<-as.numeric((freq_topic$Total_PDL1pos / freq_topic$Total))
freq_topic$Total_PDL1neg <- freq_topic$Total - freq_topic$Total_PDL1pos
freq_topic$FC <- log2(freq_topic$Total_PDL1pos/freq_topic$Total_PDL1neg)


##...Plot the frequency and the FC as bar plots...##


data_summary <- freq_topic %>%
  group_by(topic) %>%
  summarize(mean = mean(percentPDL1pos), 
            n = n(), 
            sd = sd(percentPDL1pos), 
            se = sd/sqrt(n))

line<-data_summary[data_summary$topic=="total_myeloid",]$mean

# reoder to have total myeloid first
order<-c('total_myeloid',0,1,2,3,4,5,6,7)
freq_topic$topic<-factor(freq_topic$topic, levels = order)
freq_topic<-freq_topic[order(freq_topic$topic),]


ggplot(data = freq_topic, aes(x = topic, y = percentPDL1pos, fill = topic)) + 
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3) +
  geom_hline(yintercept=line, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = '% PDL1+') +
  scale_fill_manual(values = c('#D3D3D3', colors))


# compare_means(percentPDL1pos~topic,freq_topic, method = 'wilcox.test')
# 
# 1 percentPDL1pos total_myeloid 0      0.00120    0.031   0.00120  **       Wilcoxon
# 2 percentPDL1pos total_myeloid 1      0.00562    0.12    0.00562  **       Wilcoxon
# 3 percentPDL1pos total_myeloid 2      0.436      1       0.43613  ns       Wilcoxon
# 4 percentPDL1pos total_myeloid 3      0.104      1       0.10422  ns       Wilcoxon
# 5 percentPDL1pos total_myeloid 4      0.270      1       0.27013  ns       Wilcoxon
# 6 percentPDL1pos total_myeloid 5      0.0698     0.91    0.06982  ns       Wilcoxon
# 7 percentPDL1pos total_myeloid 6      0.00000562 0.0002  5.6e-06  ****     Wilcoxon
# 8 percentPDL1pos total_myeloid 7      0.000331   0.00960 0.00033  ***      Wilcoxon


##...Produce heatmap of the percent PDL1+ for a selected population across topics...##

#by cell type
freq_topic_cell<-as.data.frame(table(topic_data$topic, topic_data$cell_type))
freq_marker_cell<-as.data.frame(table(data_marker$topic, data_marker$cell_type))
freq_topic_cell$PD.L1pos <- freq_marker_cell$Freq
names(freq_topic_cell)<-c("topic","cell_type","Total","Total_PDL1pos")
freq_topic_cell$percentPDL1pos<-as.numeric(format((freq_topic_cell$Total_PDL1pos / freq_topic_cell$Total)),digits=3)

#total myeloid
freq_topic_pooled<-as.data.frame(table(topic_data$topic))
freq_marker_pooled<-as.data.frame(table(data_marker$topic))
freq_topic_pooled$PD.L1pos <- freq_marker_pooled$Freq
names(freq_topic_pooled)<-c("topic","Total","Total_PDL1pos")
freq_topic_pooled$percentPDL1pos<-as.numeric(format((freq_topic_pooled$Total_PDL1pos / freq_topic_pooled$Total)),digits=3)

#add total myeloid data to the 
freq_topic_pooled<-add_column(freq_topic_pooled, d = 'myeloid', .after = "topic")
colnames(freq_topic_pooled)<-c("topic","cell_type","Total","Total_PDL1pos","percentPDL1pos")
freq_topic_cell<-rbind(freq_topic_pooled,freq_topic_cell)

# Plot the frequency of PDL1+ cell across all 
plot_data<-droplevels(freq_topic_cell[freq_topic_cell$cell_type=='CD163_Mac',])
plot_data<-droplevels(plot_data[plot_data$topic %in% c(1,2,3,4),])

ggplot(data = plot_data, aes(x = topic, y = percentPDL1pos, fill = topic)) + 
  geom_bar(stat='identity') +
  geom_hline(yintercept=0.7485876, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = '% PDL1+') +
  scale_fill_manual(values = c('#CA6627', '#7470AE', '#D53E88', '#74A439'))


#turn to heatmap format
cell_topic_PDL1_hmap<-dcast(freq_topic_cell, topic ~ cell_type, value.var = "percentPDL1pos")
cell_topic_PDL1_hmap<-as.matrix(cell_topic_PDL1_hmap[,-1])
cell_topic_PDL1_hmap[is.na(cell_topic_PDL1_hmap)] <- 0
rownames(cell_topic_PDL1_hmap)<-c(0,1,2,3,4,5,6,7)

heatmap.2(t(cell_topic_PDL1_hmap), 
          Colv = F, Rowv = F,
          dendrogram = 'none',
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          sepcolor="grey35",
          colsep=0:ncol(cell_topic_PDL1_hmap),
          rowsep=0:nrow(cell_topic_PDL1_hmap),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          ColSideColors = colors,
          cexRow = 1, cexCol = 2, margins = c(8,14))

# 
# # Find a CD163 example
# CD163_data<-droplevels(topic_data[topic_data$cell_type=='CD163_Mac',])
# CD163_data<-CD163_data[CD163_data$topic %in% c(1,2,3),] 
# CD163_data<-CD163_data[,c(1,2,18,63)]
# CD163_pos<-CD163_data[CD163_data$PD.L1>0.25,]
# CD163_stats<-as.data.frame(table(CD163_data$SampleID, CD163_data$topic))
# CD163_pos_stats<-as.data.frame(table(CD163_pos$SampleID, CD163_pos$topic))
# CD163_stats$PDL1pos<-CD163_pos_stats$Freq
# colnames(CD163_stats)<-c('SampleID','topic','total','total_PDL1pos')
