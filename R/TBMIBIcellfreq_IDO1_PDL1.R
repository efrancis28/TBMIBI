# TBMIBIcellfreq_IDO1_PDL1.R
# Author: Erin McCaffrey 
# Date created: 200306
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it plots the frequency of cell types across all granulomas with stacked barplots to assess the identity of PDL1
# or IDO1 positive cells across samples


library(dplyr)
library(ggplot2)
library(forcats)
library(reshape2)
library(tidyr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")

##..Drop MAC and sarcoid..##
MAC_IDs<-c(57,58,19,87)
sarcoid_IDs<-c(67,68,69,70,71,72,73,74,75,76)

data<-droplevels(data[!data$SampleID %in% c(MAC_IDs,sarcoid_IDs),])

##..Set thresholds for IDO1 and PDL1..##

IDO_thresh = 0.26
PDL1_thresh = 0.25

##..Create data frame with SampleID, cell_type, and count of that cell type in the given ROI..##

#IDO1
data_IDO1<-data[data$IDO>IDO_thresh,]
Freqs_IDO1 <- as.data.frame(table(data_IDO1$SampleID, data_IDO1$cell_type))
names(Freqs_IDO1) <- c("SampleID","cell_type","count")

#PDL1
data_PDL1<-data[data$PD.L1>PDL1_thresh,]
Freqs_PDL1<- as.data.frame(table(data_PDL1$SampleID, data_PDL1$cell_type))
names(Freqs_PDL1) <- c("SampleID","cell_type","count")

##..Get overall totals on a per sample basis..##

totals_IDO1<-aggregate(Freqs_IDO1$count, by=list(Category=Freqs_IDO1$SampleID), FUN=sum)
totals_PDL1<-aggregate(Freqs_PDL1$count, by=list(Category=Freqs_PDL1$SampleID), FUN=sum)


##..Determine frequecy of each cell type in each sample..##

#IDO1
for(i in unique(Freqs_IDO1$SampleID)) {
  frequencies_IDO1 <- Freqs_IDO1[Freqs_IDO1$SampleID==i,"count"] / totals_IDO1[totals_IDO1$Category==i,2]
  Freqs_IDO1[Freqs_IDO1$SampleID==i,"frequency"] <- frequencies_IDO1
}

#PDL1
for(i in unique(Freqs_PDL1$SampleID)) {
  frequencies_PDL1 <- Freqs_PDL1[Freqs_PDL1$SampleID==i,"count"] / totals_PDL1[totals_PDL1$Category==i,2]
  Freqs_PDL1[Freqs_PDL1$SampleID==i,"frequency"] <- frequencies_PDL1
}

##..Create stacked bars..##

cell_order<-c("CD8_T","CD4_T","CD14_Mono","imm_other","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","B_cell","CD16_CD14_Mono",
             "CD206_Mac","neutrophil","mast","Treg","CD163_Mac","gdT_cell","CD209_DC","giant_cell", "fibroblast",
             "epithelial","endothelial")

point_order_gran<-c(64, 65, 21, 84, 42, 88, 28, 89, 85, 13, 35, 36, 14, 15, 
                    6, 7, 33, 34, 26, 27, 40, 61, 47, 48, 54, 55)


##..Create stacked bar plot frequency across ROIs..##


#PDL1
Freqs_PDL1$cell_type<- factor(Freqs_PDL1$cell_type, levels=rev(cell_order))
Freqs_PDL1<-Freqs_PDL1[order(Freqs_PDL1$cell_type),]
Freqs_PDL1$SampleID <- factor(Freqs_PDL1$SampleID, levels=point_order_gran)
Freqs_PDL1<-Freqs_PDL1[order(Freqs_PDL1$SampleID),]

colorkey<-read.csv('colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% cell_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = cell_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

PDL1_bar<-ggplot(Freqs_PDL1, aes(x=SampleID, y=frequency, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values =rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  theme(legend.position = 'none') +
  guides(fill=guide_legend(title="Cell Type"))
PDL1_bar


PDL1_count_bar<-ggplot(Freqs_PDL1, aes(x=SampleID, y=count, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values =rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  theme(legend.position = 'none') +
  guides(fill=guide_legend(title="Cell Type"))
PDL1_count_bar


#IDO1
Freqs_IDO1$cell_type<- factor(Freqs_IDO1$cell_type, levels=rev(cell_order))
Freqs_IDO1<-Freqs_IDO1[order(Freqs_IDO1$cell_type),]
Freqs_IDO1$SampleID <- factor(Freqs_IDO1$SampleID, levels=point_order_gran)
Freqs_IDO1<-Freqs_IDO1[order(Freqs_IDO1$SampleID),]


IDO1_bar<-ggplot(Freqs_IDO1, aes(x=SampleID, y=frequency, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values =rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  theme(legend.position = 'none') +
  guides(fill=guide_legend(title="Cell Type"))
IDO1_bar


IDO1_count_bar<-ggplot(Freqs_IDO1, aes(x=SampleID, y=count, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values =rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  theme(legend.position = 'none') +
  guides(fill=guide_legend(title="Cell Type"))
IDO1_count_bar

