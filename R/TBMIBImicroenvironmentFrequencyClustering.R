# TBMIBImicroenvironmentFrequencyClustering.R
# Author: Erin McCaffrey 
# Date created: 200520
# Overview: This script creates a matrix of all patients and the frequency of each ME in that patient. 
# It next runs a clustering analysis of that data based on correlation and euclidean distance with 
# heirarchical clustering.

library(reshape2)
library(corrplot)
library(psych)
library(gplots)
library(ggplot2)
library(colorspace)
library("GMD")

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##...Enumerate frequency of each topic per patient..##

topic_freqs<-as.data.frame(table(topic_data$SampleID, topic_data$topic))
colnames(topic_freqs)<-c('SampleID','topic','count')
patient_totals<-as.vector(table(topic_data$SampleID))
topic_freqs$total<-rep(patient_totals, n=8)
topic_freqs$freq<-topic_freqs$count / topic_freqs$total

##..Cast to a matrix that is all patients (columns) and all topics (rows)..##

topic_mat<-reshape2::dcast(topic_freqs, topic ~ SampleID, value.var = "freq")
topicmat.m<-as.matrix(topic_mat[,-1])
row.names(topicmat.m)<-topic_mat$topic

##..Plot heatmap with hierarchical clustering along columns based on euclidean..##

colors<-c('#4B9B79', '#CA6627', '#7470AE', '#D53E88', '#74A439', '#DDAE3B', '#BB2D34', '#668AB7')

heatmap.2(topicmat.m, 
          Colv = T, Rowv = F,
          hclustfun = hclust,
          scale = "none",
          dendrogram = "column",
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(topicmat.m),
          rowsep=0:nrow(topicmat.m),
          RowSideColors = colors,
          breaks=seq(0, 0.5, length.out=101),
          cexRow = 2, cexCol = 1, margins = c(8,14))

##..Cluster analysis of hclust with euclidean..##

# get ideal number of clusters with wss
data<-topicmat.m
dist_mat<-dist(t(data), method = 'euclidean')
hc<-hclust(dist(t(data), method = 'euclidean'), method = 'complete')
plot(hc)
rect.hclust(hc, k = 5)

cluster_stats<-css.hclust(dist_mat, hclust.obj=hc, k=10)
plot(attr(cluster_stats,"meta")$hclust.obj)

# plot explained variance

ggplot(cluster_stats, aes(x=k, y=ev)) +
  geom_line() +
  geom_point()

##..Plot heatmap with hierarchical clustering along columns based on correlation..##

heatmap.2(topicmat.m, 
          Colv = T, Rowv = F,
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="complete"),
          scale = "none",
          dendrogram = "column",
          trace = "none",
          col = rev(sequential_hcl(100, palette = 'Reds 2')),
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(topicmat.m),
          rowsep=0:nrow(topicmat.m),
          RowSideColors = colors,
          breaks=seq(0, 1, length.out=101),
          cexRow = 2, cexCol = 1, margins = c(8,14))

##..Cluster analysis of Pearson correlation-based clustering..##

data<-topicmat.m
corr<-corr.test(data, y = NULL, use = "pairwise",method="pearson",adjust="fdr")

# get ideal number of clusters with wss
data_corr<-corr$r
dist_mat<-as.dist(1 - data_corr)
hc <- hclust(dist_mat, method = 'complete')
plot(hc)
rect.hclust(hc, k = 5)

cluster_stats<-css.hclust(dist_mat, hclust.obj=hc, k=10)
plot(attr(cluster_stats,"meta")$hclust.obj)


# plot explained variance

ggplot(cluster_stats, aes(x=k, y=ev)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0,10,1)) 


