# TBMIBIquantifyIDOPosFreqAcrossMEsAndCells.R
# Author: Erin McCaffrey
# Date created: 200106

library(ggpubr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(forcats)
library(tibble)
library(reshape2)
library(gplots)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Filter non-myeloid cell types...##

cell_types<-c('myeloid')
topic_data <- droplevels(topic_data[topic_data$cell_lin %in% cell_types, ])


##...Get percent of cells positive for IDO and ratio of pos:neg...##

marker_thresh <- 0.26
data_marker<-droplevels(topic_data[topic_data$IDO>marker_thresh, ])

#Across samples
freq_sample<-as.data.frame(table(topic_data$SampleID, topic_data$topic, topic_data$cell_type))
freq_marker_sample<-as.data.frame(table(data_marker$SampleID, data_marker$topic, data_marker$cell_type))
freq_sample$IDOpos <- freq_marker_sample$Freq
names(freq_sample)<-c("SampleID","topic","cell_type","Total","Total_IDOpos")
freq_sample$Total_IDOneg <- freq_sample$Total - freq_sample$Total_IDOpos
freq_sample$percentIDOpos<-as.numeric(format((freq_sample$Total_IDOpos / freq_sample$Total)*100),digits=3)
freq_sample$FC <- log2(freq_sample$Total_IDOpos/freq_sample$Total_IDOneg)

##...Visualize plot of FC across topics (all cell types, boxplot...##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

ggplot(data = freq_sample, aes(x = topic, y = FC, fill = topic)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = 'log2(IDO+ / IDO-)') +
  scale_fill_manual(values = colors)

##...Produce table of frequency of IDO+ for all myeloid cells and then each topic...##

# across topics and points
freq_topic<-as.data.frame(table(topic_data$topic, topic_data$SampleID))
freq_marker_topic<-as.data.frame(table(data_marker$topic, data_marker$SampleID))
freq_topic$IDOpos <- freq_marker_topic$Freq
names(freq_topic)<-c("topic","SampleID","Total","Total_IDOpos")

# total myeloid
freq_myeloid<-as.data.frame(table(topic_data$SampleID))
freq_myeloid<-add_column(freq_myeloid, topic = 'total_myeloid', .before = "Var1")
freq_myeloid_marker<-as.data.frame(table(data_marker$SampleID))
freq_myeloid$IDOpos <- freq_myeloid_marker$Freq
names(freq_myeloid)<-c("topic","SampleID","Total","Total_IDOpos")

# merge and get frequency
freq_topic<-rbind(freq_myeloid, freq_topic)
freq_topic$percentIDOpos<-as.numeric((freq_topic$Total_IDOpos / freq_topic$Total))
freq_topic$Total_IDOneg <- freq_topic$Total - freq_topic$Total_IDOpos
freq_topic$FC <- log2(freq_topic$Total_IDOpos/freq_topic$Total_IDOneg)


##...Plot the frequency and the FC as bar plots...##


data_summary <- freq_topic %>%
  group_by(topic) %>%
  summarize(mean = mean(percentIDOpos), 
            n = n(), 
            sd = sd(percentIDOpos), 
            se = sd/sqrt(n))

line<-data_summary[data_summary$topic=="total_myeloid",]$mean

# reoder to have total myeloid first
order<-c('total_myeloid',0,1,2,3,4,5,6,7)
freq_topic$topic<-factor(freq_topic$topic, levels = order)
freq_topic<-freq_topic[order(freq_topic$topic),]

ggplot(data = freq_topic, aes(x = topic, y = percentIDOpos, fill = topic)) + 
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3) +
  geom_hline(yintercept=line, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = '% IDO+') +
  scale_fill_manual(values = c('#D3D3D3', colors))

# compare_means(percentIDOpos~topic,plot_data, method = 'wilcox.test')
# 
# 1 percentIDOpos total_myeloid 0      0.0930   1    0.09295  ns       Wilcoxon
# 2 percentIDOpos total_myeloid 1      0.0276   0.74 0.02757  *        Wilcoxon
# 3 percentIDOpos total_myeloid 2      0.745    1    0.74531  ns       Wilcoxon
# 4 percentIDOpos total_myeloid 3      0.541    1    0.54071  ns       Wilcoxon
# 5 percentIDOpos total_myeloid 4      0.295    1    0.29512  ns       Wilcoxon
# 6 percentIDOpos total_myeloid 5      0.218    1    0.21788  ns       Wilcoxon
# 7 percentIDOpos total_myeloid 6      0.00393  0.13 0.00393  **       Wilcoxon
# 8 percentIDOpos total_myeloid 7      0.0139   0.42 0.01391  *        Wilcoxon


##...Produce heatmap of the percent IDO+ for each population across topics...##

#by cell type
freq_topic_cell<-as.data.frame(table(topic_data$topic, topic_data$cell_type))
freq_marker_cell<-as.data.frame(table(data_marker$topic, data_marker$cell_type))
freq_topic_cell$IDOpos <- freq_marker_cell$Freq
names(freq_topic_cell)<-c("topic","cell_type","Total","Total_IDO1pos")
freq_topic_cell$percentIDO1pos<-as.numeric(format((freq_topic_cell$Total_IDO1pos / freq_topic_cell$Total)),digits=3)

#total myeloid
freq_topic_pooled<-as.data.frame(table(topic_data$topic))
freq_marker_pooled<-as.data.frame(table(data_marker$topic))
freq_topic_pooled$IDOpos <- freq_marker_pooled$Freq
names(freq_topic_pooled)<-c("topic","Total","Total_IDO1pos")
freq_topic_pooled$percentIDO1pos<-as.numeric(format((freq_topic_pooled$Total_IDO1pos / freq_topic_pooled$Total)),digits=3)

#add total myeloid data to the 
freq_topic_pooled<-add_column(freq_topic_pooled, d = 'myeloid', .after = "topic")
colnames(freq_topic_pooled)<-c("topic","cell_type","Total","Total_IDO1pos","percentIDO1pos")
freq_topic_cell<-rbind(freq_topic_pooled,freq_topic_cell)

# Plot the frequency of PDL1+ cell across all CD11b/c
plot_data<-droplevels(freq_topic_cell[freq_topic_cell$cell_type=='CD11c_DC/Mono',])
plot_data<-droplevels(plot_data[plot_data$topic %in% c(1,2,3,4,5),])

ggplot(data = plot_data, aes(x = topic, y = percentIDO1pos, fill = topic)) + 
  geom_bar(stat='identity') +
  geom_hline(yintercept=0.5345599, linetype="dashed", color = "black", size =1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x = 'Topic') + 
  labs(y = '% IDO1+') +
  scale_fill_manual(values = c('#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B'))

#turn to heatmap format
cell_topic_IDO1_hmap<-dcast(freq_topic_cell, topic ~ cell_type, value.var = "percentIDO1pos")
cell_topic_IDO1_hmap<-as.matrix(cell_topic_IDO1_hmap[,-1])
cell_topic_IDO1_hmap[is.na(cell_topic_IDO1_hmap)] <- 0
rownames(cell_topic_IDO1_hmap)<-c(0,1,2,3,4,5,6,7)

heatmap.2(t(cell_topic_IDO1_hmap), 
          Colv = F, Rowv = F,
          dendrogram = 'none',
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          sepcolor="grey35",
          colsep=0:ncol(cell_topic_IDO1_hmap),
          rowsep=0:nrow(cell_topic_IDO1_hmap),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          ColSideColors = colors,
          cexRow = 1, cexCol = 2, margins = c(8,14))


# Find a good CD163 example
CD11c_data<-droplevels(topic_data[topic_data$cell_type=='CD11c_DC/Mono',])
CD11c_data<-CD11c_data[CD11c_data$topic %in% c(1,2,3,4,5),] 
CD11c_data<-CD11c_data[,c(1,2,29,63)]
CD11c_pos<-CD11c_data[CD11c_data$IDO>0.26,]
CD11c_stats<-as.data.frame(table(CD11c_data$SampleID, CD11c_data$topic))
CD11c_pos_stats<-as.data.frame(table(CD11c_pos$SampleID, CD11c_pos$topic))
CD11c_stats<-droplevels(CD11c_stats[CD11c_stats$Var1 %in% unique(CD11c_pos_stats$Var1),])
CD11c_stats$PDL1pos<-CD11c_pos_stats$Freq
colnames(CD11c_stats)<-c('SampleID','topic','total','total_IDO1pos')
