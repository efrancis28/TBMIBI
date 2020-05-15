# TBMIBIcellfreq.R
# Author: Erin McCaffrey 
# Date created: 190318
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it plots the frequency of cell types across all granulomas with boxplots and stacked barplots. Done for both major
# lineage and then within the immune cell compartment.


library(dplyr)
library(viridis)
library(ggplot2)
library(forcats)
library(reshape2)
library(tidyr)

##..Import data..##
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")

data_lineage<-data %>% select(SampleID,lineage)

data_immlineage<-data %>% select(SampleID,cell_lin)
data_immlineage<-droplevels(data_immlineage[!data_immlineage$cell_lin=="nonimmune",])

data_immune<-data %>% select(SampleID,cell_type)
data_immune<-droplevels(data_immune[!data_immune$cell_type  %in% c('epithelial','endothelial','fibroblast'),])


##..Create data frame with SampleID, cell_type, and count of that cell type in the given ROI..##

#lineage
Freqs_lineage <- table(data_lineage$SampleID, data_lineage$lineage)
Freqs_lineage <- as.data.frame(Freqs_lineage)
names(Freqs_lineage) <- c("SampleID","lineage","count")

#immune
Freqs_immune <- table(data_immune$SampleID, data_immune$cell_type)
Freqs_immune <- as.data.frame(Freqs_immune)
names(Freqs_immune) <- c("SampleID","cell_type","count")

#all cell types (total frequency)
Freqs_total <- table(data$SampleID, data$cell_type)

Freqs_total <- as.data.frame(Freqs_total)
names(Freqs_total) <- c("SampleID","cell_type","count")

#immune lineage
Freqs_immunelin <- table(data_immlineage$SampleID, data_immlineage$cell_lin)
Freqs_immunelin <- as.data.frame(Freqs_immunelin)
names(Freqs_immunelin) <- c("SampleID","cell_lin","count")

##..Get overall totals on a per sample basis..##
totals_lineage<-aggregate(Freqs_lineage$count, by=list(Category=Freqs_lineage$SampleID), FUN=sum)
totals_immune<-aggregate(Freqs_immune$count, by=list(Category=Freqs_immune$SampleID), FUN=sum)
totals_all<-aggregate(Freqs_total$count, by=list(Category=Freqs_total$SampleID), FUN=sum)
total_immune_lineage<-aggregate(Freqs_immunelin$count, by=list(Category=Freqs_immunelin$SampleID), FUN=sum)

##..Drop MAC..##
MAC_IDs<-c(57,58,19,87)

##..Determine frequecy of each cell type in each sample..##

#lineage
for(i in unique(Freqs_lineage$SampleID)) {
  frequencies <- Freqs_lineage[Freqs_lineage$SampleID==i,"count"] / totals_lineage[totals_lineage$Category==i,2]
  Freqs_lineage[Freqs_lineage$SampleID==i,"frequency"] <- frequencies
}

#immune
for(i in unique(Freqs_immune$SampleID)) {
  frequencies_imm <- Freqs_immune[Freqs_immune$SampleID==i,"count"] / totals_immune[totals_immune$Category==i,2]
  Freqs_immune[Freqs_immune$SampleID==i,"frequency"] <- frequencies_imm
}

#all cell types
for(i in unique(Freqs_total$SampleID)) {
  frequencies_total <- Freqs_total[Freqs_total$SampleID==i,"count"] / totals_all[totals_all$Category==i,2]
  Freqs_total[Freqs_total$SampleID==i,"frequency"] <- frequencies_total
}

#immune lineage
for(i in unique(Freqs_immunelin$SampleID)) {
  frequencies_immlin <- Freqs_immunelin[Freqs_immunelin$SampleID==i,"count"] / total_immune_lineage[total_immune_lineage$Category==i,2]
  Freqs_immunelin[Freqs_immunelin$SampleID==i,"frequency"] <- frequencies_immlin
}


##..Create boxplot across cell types..##

imm_order<-c("CD8_T","CD4_T","CD14_Mono","imm_other","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","B_cell","CD16_CD14_Mono",
             "CD206_Mac","neutrophil","mast","Treg","CD163_Mac","gdT_cell","CD209_DC","giant_cell")
lin_order<-c("immune","endothelial","fibroblast","epithelial")
immlin_order<-c("myeloid","lymphocyte","other","granulocyte")
point_order_gran<-c(64, 65, 21, 84, 42, 88, 28, 89, 85, 13, 35, 36, 14, 15, 57, 58, 19, 87, 6, 7, 33, 34, 26, 27, 40, 61, 47, 48, 54, 55)
point_order_sarc<-c(67,68,69,70,71,72,73,74,75,76)
point_order_all<-c(point_order_gran,point_order_sarc)
x_names <- c("1-1","1-2","2-1","2-2","3-1","3-2","4-1","4-2","5-1","5-2","6-1","6-2","7-1","7-2","8-1","8-2","9-1","9-2","10-1","10-2","11-1","11-2","12-1","12-2","13-1","13-2","14-1","14-2","15-1","15-2")

##..Create a barplot of cell totals across FOVs..##

plot_data<-droplevels(totals_all[!totals_all$Category %in% MAC_IDs & !totals_all$Category %in% point_order_sarc,])

cell_bar<-ggplot(plot_data, aes(x=reorder(as.factor(Category),-x), y=x)) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_manual("grey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) 
cell_bar

#reorder by cell type and point order
Freqs_lineage$SampleID <- factor(Freqs_lineage$SampleID, levels=point_order_all)
Freqs_immune$SampleID <- factor(Freqs_immune$SampleID, levels=point_order_all)
Freqs_total$SampleID <- factor(Freqs_total$SampleID, levels=point_order_all)
Freqs_immunelin$SampleID <- factor(Freqs_immunelin$SampleID, levels=point_order_all)

Freqs_lineage$lineage <- factor(Freqs_lineage$lineage, levels=lin_order)
Freqs_immune$cell_type <- factor(Freqs_immune$cell_type, levels=rev(imm_order))
Freqs_immunelin$cell_lin <- factor(Freqs_immunelin$cell_lin, levels=immlin_order)

#separate the mycobacteria and sarcoid data from each other
Freqs_lineage_gran<-Freqs_lineage[Freqs_lineage$SampleID %in% point_order_gran,] #lineage mycobacteria
Freqs_lineage_sarc<-Freqs_lineage[Freqs_lineage$SampleID %in% point_order_sarc,] #lineage sarcoid

Freqs_immune_gran<-Freqs_immune[Freqs_immune$SampleID %in% point_order_gran,] #immune mycobacteria
Freqs_immune_sarc<-Freqs_immune[Freqs_immune$SampleID %in% point_order_sarc,] #immune sarcoid

Freqs_total_gran<-Freqs_total[Freqs_total$SampleID %in% point_order_gran,] #all mycobacteria
Freqs_total_sarc<-Freqs_total[Freqs_total$SampleID %in% point_order_sarc,] #all sarcoid

Freqs_immunelin_gran<-Freqs_immunelin[Freqs_immunelin$SampleID %in% point_order_gran,] #immune mycobacteria
Freqs_immunelin_sarc<-Freqs_immunelin[Freqs_immunelin$SampleID %in% point_order_sarc,] #immune sarcoid

##..Create box plots of pooled ROIs..##

lin_box<-ggplot(Freqs_lineage_gran[!Freqs_lineage_gran$SampleID %in% MAC_IDs,], aes(x=lineage, y=frequency, fill=lineage)) + geom_boxplot() +
  theme_bw() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
lin_box

Freqs_immune_sarc$cell_type<-reorder(Freqs_immune_sarc$cell_type, new.order=imm_order)

colorkey<-read.csv('colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% imm_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = imm_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color<-as.vector(colorkey_imm$code)

imm_box<-ggplot(Freqs_immune_gran[!Freqs_immune_gran$SampleID %in% MAC_IDs,], aes(x=fct_reorder(cell_type,frequency,.fun=median,.desc=TRUE), y=frequency, fill=cell_type)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank(), axis.text.x = element_text(angle=35,hjust=1)) +
  labs(x="Cell Type") +
  labs(y="Frequency") +
  ggtitle("Frequency of Immune Cell Types") +
  guides(fill=guide_legend(title="Cell Type"))
imm_box

##..Create bar plots lineage and immune for myco of pooled ROIs..##

lin_bar<-ggplot(Freqs_lineage_gran[!Freqs_lineage_gran$SampleID %in% MAC_IDs,], aes(x=lineage, y=frequency)) +
  theme_bw() +
  scale_fill_manual("white") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) +
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3)
lin_bar

immlin_bar<-ggplot(Freqs_immunelin_gran[!Freqs_immunelin_gran$SampleID %in% MAC_IDs,], aes(x=cell_lin, y=frequency)) +
  theme_bw() +
  scale_fill_manual("white") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.ticks.x=element_blank()) +
  stat_summary(geom = "bar", fun.y = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.3)
immlin_bar

##..Create stacked bar plot across ROIs..##

Freqs_lineage_gran$lineage <- factor(Freqs_lineage_gran$lineage, levels=rev(lin_order))
Freqs_lineage_gran$SampleID <- factor(Freqs_lineage_gran$SampleID, levels=point_order_gran)
Freqs_lineage_gran<-Freqs_lineage_gran[order(Freqs_lineage_gran$SampleID),]

Freqs_lineage_sarc$lineage <- factor(Freqs_lineage_sarc$lineage, levels=rev(lin_order))

lin_bar<-ggplot(Freqs_lineage_gran[!Freqs_lineage_gran$SampleID %in% MAC_IDs,], aes(x=SampleID, y=frequency, fill=lineage)) + 
  theme_bw() +
  scale_fill_manual(values =c('#FF9966','#99FF99','#CC3333','#3366FF')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_discrete(breaks=point_order_gran,labels=x_names) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
lin_bar


#re-order to have the samples in order specified and immune phenos in order of decreasing median
Freqs_immune_gran$SampleID <- factor(Freqs_immune_gran$SampleID, levels=point_order_gran)
Freqs_immune_gran<-Freqs_immune_gran[order(Freqs_immune_gran$SampleID),]

Freqs_immune_gran$cell_type <- factor(Freqs_immune_gran$cell_type, levels=rev(imm_order))
Freqs_immune_gran<-Freqs_immune_gran[order(Freqs_immune_gran$cell_type),]

Freqs_immune_sarc$cell_type <- factor(Freqs_immune_sarc$cell_type, levels=rev(imm_order))
Freqs_immune_sarc<-Freqs_immune_sarc[order(Freqs_immune_sarc$cell_type),]

#plot just the TB samples with heirarchical order based on the pearson correlation 
Freqs_immune_TB<-droplevels(Freqs_immune_gran[!Freqs_immune_gran$SampleID %in% MAC_IDs,])
h_order<-c(33,34,47,48,6,7,61,26,27,15,21,84,40,14,88,42,13,89,85,64,65,55,35,54,36,28)
Freqs_immune_TB$SampleID <- factor(Freqs_immune_TB$SampleID, levels=h_order)
Freqs_immune_TB<-Freqs_immune_TB[order(Freqs_immune_TB$SampleID),]

imm_bar<-ggplot(Freqs_immune_TB, aes(x=SampleID, y=frequency, fill=cell_type)) + 
  theme_bw() +
  scale_fill_manual(values = rev(color)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #scale_x_discrete(breaks=h_order,labels=h_order) +
  geom_bar(stat = "identity", colour="white", lwd=0.2) +
  xlab("Sample ID (Patient-ROI)") +
  ylab("Frequency") +
  guides(fill=guide_legend(title="Cell Type"))
imm_bar

##..Append sample type info and save frequency data..##

Freqs_immune<-Freqs_immune %>% mutate(Tissue=case_when(Freqs_immune$SampleID %in% point_order_gran~ "mycobacteria",
                                                       Freqs_immune$SampleID %in% point_order_sarc ~ "sarcoid"))


Freqs_lineage<-Freqs_lineage %>% mutate(Tissue=case_when(Freqs_lineage$SampleID %in% point_order_gran~ "mycobacteria",
                                                       Freqs_lineage$SampleID %in% point_order_sarc ~ "sarcoid"))

Freqs_total<-Freqs_total %>% mutate(Tissue=case_when(Freqs_total$SampleID %in% point_order_gran~ "mycobacteria",
                                                     Freqs_total$SampleID %in% point_order_sarc ~ "sarcoid"))

Freqs_immunelin<-Freqs_immunelin %>% mutate(Tissue=case_when(Freqs_immunelin$SampleID %in% point_order_gran~ "mycobacteria",
                                                             Freqs_immunelin$SampleID %in% point_order_sarc ~ "sarcoid"))


write.csv(Freqs_immune, file="immune_cell_freqs.csv",row.names = FALSE)
write.csv(Freqs_lineage, file="lineage_freqs.csv",row.names = FALSE)
write.csv(Freqs_total, file="total_freqs.csv",row.names = FALSE)
write.csv(Freqs_immunelin, file="immlineage_freqs.csv",row.names = FALSE)

write.csv(totals_lineage, file="lineage_cell_totals.csv",row.names = FALSE)
write.csv(totals_immune, file="immune_cell_totals.csv",row.names = FALSE)


# ##..Custom color scale..##
# 
# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# 
# color<-gg_color_hue(17)
# colorkey<-data.frame(imm_order)
# colorkey$code<-color
# write.csv(colorkey,file='colorkey_R.csv',row.names = FALSE)

