# TBMIBIquantifyExpressionAcrossAllMEs.R
# Author: Erin McCaffrey
# Date created: 191206
# Overview: This script reads in the csv for the topic loadings for each cell. It produces
# an expression heatmap based on the assignment of each cell to a topic (composite) and by multiplying
# each cell's expression by the weights for all topics. 

library(plyr)
library(dplyr)
library(pals)
library(colorspace)
library(gplots)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##..Define markers and topics..##

major<-c('CD45','SMA','E.cadherin','CD31','Keratin.pan',"Vimentin")
immune<-c("CD3","CD20","CD14","CD16","CD11c","HLA.DR.DQ.DP","MastChyTry","MPO")
myeloid<-c("CD14","CD16","CD11c","CD11b","HLA.DR.DQ.DP","CD68","CD163","CD206","CD209")
lymph<-c("CD3","CD20","CD4","CD8","Foxp3")
functional<-c("IDO","PD.L1","iNOS","CD36", "H3K9Ac", "PD.1", "Lag3", "HLA.Class.1", "Collagen.1",'IFNg', 'CD103','Ki67')

heatmap_channels<-unique(c(major,immune,myeloid,lymph, functional))
topics<-c(0,1,2,3,4,5,6,7)

##..Produce heatmap matrix based on assignment to topics..##

hm_allclusters <- matrix(, nrow = length(topics), ncol = length(heatmap_channels))
for(i in 1:length(topics)) {
  temp_mat <- topic_data[topic_data$topic==topics[i], heatmap_channels]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}


rownames(hm_allclusters) <- paste("topic_", topics, sep = "")
colnames(hm_allclusters) <- heatmap_channels
hm_allclusters

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')
heatmap.2(hm_allclusters, 
          Colv = T, Rowv = T,
          hclustfun = hclust,
          scale = "column",
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = as.vector(ocean.curl(100)),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          RowSideColors = colors,
          cexRow = 2, cexCol = 2, margins = c(8,14))

##..Produce heatmap matrix without assignments of cells to topics (weighted expression)..##

hm_allclusters_prob <- matrix(, nrow = length(topics), ncol = length(heatmap_channels))
exprs_mat <- topic_data[, heatmap_channels] 
for(i in 1:length(topics)) {
  topic_weight <-matrix(topic_data[,paste("Topic.", topics[i], sep = "")])
  temp_mat <- sapply(exprs_mat, '*', topic_weight)
  hm_allclusters_prob[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

rownames(hm_allclusters_prob) <- paste("topic_", topics, sep = "")
colnames(hm_allclusters_prob) <- heatmap_channels
hm_allclusters_prob

heatmap.2(hm_allclusters_prob[,23:34], 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "column",
          dendrogram = "column",
          trace = "none",
          col = diverging_hcl(100, palette = 'Blue-Red 3'),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters),
          rowsep=0:nrow(hm_allclusters),
          RowSideColors = colors,
          cexRow = 2, cexCol = 2, margins = c(8,14))

