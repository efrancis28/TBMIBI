# TBMIBIquantifyCellTypeLoadingsAcrossMEs.R
# Author: Erin McCaffrey
# Date created: 191211
# Overview: This script reads in the csv for the topic loadings for each cell

library(dplyr)
library(gplots)
require(pals)
library(colorspace)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Define topics and cell_types...##

phenos<-as.vector(unique(topic_data$cell_type))
topics<-as.vector(colnames(topic_data[,55:62]))

##...For each topic and cell type, produce heatmap with the mean loading...##

# Produce heatmap matrix

hm_alltopics <- matrix(, nrow = length(phenos) , ncol=length(topics))
for(i in 1:length(phenos)) {
  temp_mat <- topic_data[topic_data$cell_type==phenos[i], topics]
  hm_alltopics[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}


rownames(hm_alltopics) <- phenos
colnames(hm_alltopics) <- topics
hm_alltopics

# Visualize

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

heatmap.2(hm_alltopics, 
          Colv = F, Rowv = T,
          hclustfun = hclust,
          scale = "row",
          dendrogram = c("row"),
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          sepcolor="grey35",
          colsep=0:ncol(hm_alltopics),
          rowsep=0:nrow(hm_alltopics),
          sepwidth=c(0.01,0.01),
          symkey=F,
          density.info = 'none',
          key.title = '',
          ColSideColors = colors,
          cexRow = 1, cexCol = 1, margins = c(8,14))


