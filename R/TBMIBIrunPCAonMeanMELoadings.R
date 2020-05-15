# TBMIBIrunPCAonMeanMELoadings.R
# Author: Erin McCaffrey
# Date created: 200118
# Overview: This script reads in the csv for the topic loadings for each cell across all TB patients. 
# It determines the mean topic loading per patient, runs PCA on this data, and visualizes a PCA colored
# by patient ID.

library(ggplot2)
library(ggrepel)
library(factoextra)

##...Load in data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/topic_analysis")
topic_data<-read.csv('all_TB_topic_annotation.csv')

##..Produce matrix of all topic loadings versus all samples..##

SampleIDs<-unique(topic_data$SampleID)
topics<-c(0,1,2,3,4,5,6,7)
topics<- paste("Topic.", topics, sep = "")

##..Produce matrix of mean topic loadings..##

alltopics.mat <- matrix(, nrow = length(topics), ncol = length(SampleIDs))
for(i in 1:length(SampleIDs)) {
  temp_mat <- topic_data[topic_data$SampleID==SampleIDs[i], topics]
  alltopics.mat[,i] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}


rownames(alltopics.mat) <- topics
colnames(alltopics.mat) <- SampleIDs
alltopics.mat<-as.data.frame(t(alltopics.mat))
alltopics.mat$SampleID<-rownames(alltopics.mat)
alltopics.mat$PatientROI<-c("8-1","8-2","5-2","7-1","7-2","2-1","10-1","10-2","4-1","9-1","9-2",
                            "6-1","6-2","11-1","3-1","12-1","12-2","13-1","13-2","11-2",
                            "1-1","1-2","2-2","5-1","3-2","4-2")



##..Run and visualize PCA..##

gran.topic.pca<-prcomp(alltopics.mat[,1:8], center=T)
fviz_eig(gran.topic.pca)
alltopics.mat$PC1<-gran.topic.pca$x[,1]
alltopics.mat$PC2<-gran.topic.pca$x[,2]

ggplot(alltopics.mat, aes(PC1,PC2, label=PatientROI)) +
  # geom_point(aes(x=-0.5370, y=-0.1732), colour="red", size =6) +
  # geom_point(aes(x=0.0270, y=0.4066), colour="green", size = 6) +
  # geom_point(aes(x=0.4191, y=-0.2313), colour="blue", size=6) +
  geom_point(size=2) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black", size =0.5) +
  geom_vline(xintercept=0, linetype="dashed", color = "black", size =0.5) +
  geom_label_repel(aes(label=PatientROI)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="PC 1 (45% of Variance Explained)") +
  labs(y="PC 2 (23% of Variance Explained)") +
  labs(color = "Patient-ROI") +
  ggtitle("PCA Topics")

