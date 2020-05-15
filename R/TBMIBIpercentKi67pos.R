# TBMIBIpercentKi67pos.R
# Author: Erin McCaffrey 
# Date created: 190501
# Overview: This script reads in the csv for annotated, cell-size normalized, scaled by 1000
# expression data. It asinh transforms the data, then looks at the percent of subsets that
# are Ki67 positive. 

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(exactRankTests)

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")

data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN'),])

##..Get just the lymphocytes..##

lymphs<-c("B_cell","CD4_T","CD8_T","Treg")

##..Get the percent Ki67 pos cells in each subset..#

data_lymph<-droplevels(data_myco[data_myco$cell_type %in% lymphs, ])
data_lymph_Ki67<-droplevels(data_lymph[data_lymph$Ki67>0.27578920, ])


##..Get appropriate color key..##

colorkey<-read.csv("colorkey_R.csv")
lymph_color<-droplevels(colorkey[colorkey$imm_order %in% lymphs,])
lymph_color$imm_order<-factor(lymph_color$imm_order, levels = lymphs)
lymph_color<-lymph_color[order(lymph_color$imm_order),]
color<-as.vector(lymph_color$code)

#bulk
freq<-as.data.frame(table(data_lymph$cell_type))
freq_Ki67<-as.data.frame(table(data_lymph_Ki67$cell_type))
freq$Ki67pos<-freq_Ki67$Freq
names(freq)<-c("Cell_Type","Total","Total_Ki67pos")
freq$percent<-as.numeric(format((freq$Total_Ki67pos/freq$Total*100), digits=3))

bulk<-ggplot(data=freq, aes(x=reorder(Cell_Type, -percent), y=percent, color=Cell_Type)) +
  geom_text(aes(label=percent), vjust=-0.3, size=3.5)+
  geom_bar(stat="identity", fill="white", lwd=1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="% Ki67+ of Total") +
  guides(fill=FALSE) 
bulk

#Across samples
freq_sample<-as.data.frame(table(data_lymph$cell_type, data_lymph$SampleID))
freq_Ki67_sample<-as.data.frame(table(data_lymph_Ki67$cell_type,data_lymph_Ki67$SampleID))
freq_sample$Ki67pos<-0
freq_sample[freq_sample$Var2 %in% freq_Ki67_sample$Var2 & freq_sample$Var1 %in% freq_Ki67_sample$Var1,]$Ki67pos<-freq_Ki67_sample$Freq
names(freq_sample)<-c("Cell_Type","SampleID","Total","Total_Ki67pos")
freq_sample$percent<-as.numeric(format((freq_sample$Total_Ki67pos / freq_sample$Total)*100),digits=3)


my_comparisons <- list(c("Treg","CD4_T"),c("Treg","CD8_T"),c("Treg","B_cell"))
#plot
per_sample<-ggplot(na.omit(freq_sample), aes(x=fct_reorder(Cell_Type,percent,.fun=median,.desc=TRUE), y=percent, fill=Cell_Type)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.x=element_blank()) +
  stat_compare_means(data=na.omit(freq_sample),mapping=aes(x=fct_reorder(Cell_Type,percent,.fun=median,.desc=TRUE), y=percent)
                     ,label = "p.format", method= "wilcox.test",
                     comparisons=my_comparisons) +
  theme(legend.position = "none") +
  labs(x="Cell Type") +
  labs(y="% Ki67+ of Total") +
  guides(fill=guide_legend(title="Cell Type"))
per_sample

compare_means(percent ~ Cell_Type, freq_sample, method= "wilcox.test")

wilcox.exact(percent ~ Cell_Type, freq_sample, Cell_Type %in% c('Treg','B_cell'))

# .y.     group1 group2        p   p.adj p.format p.signif method  
# <chr>   <chr>  <chr>     <dbl>   <dbl> <chr>    <chr>    <chr>   
# 1 percent B_cell CD4_T  0.0109   0.033   0.01086  *        Wilcoxon
# 2 percent B_cell CD8_T  0.0701   0.14    0.07010  ns       Wilcoxon
# 3 percent B_cell Treg   0.000105 0.00063 0.00011  ***      Wilcoxon
# 4 percent CD4_T  CD8_T  0.253    0.25    0.25265  ns       Wilcoxon
# 5 percent CD4_T  Treg   0.000333 0.0013  0.00033  ***      Wilcoxon
# 6 percent CD8_T  Treg   0.000238 0.00120 0.00024  ***      Wilcoxon
