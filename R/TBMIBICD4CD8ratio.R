# TBMIBICD4CD8ratio.R
# Author: Erin McCaffrey 
# Date created: 190501
# Overview: This script reads in the frequency data for granuloma immune cell population. It then determines
# the ratio of CD4 T cells to CD8 T cells (log2(CD4/CD8)) and plots them in descending order (ie. more CD4s on 
# the left side and less CD4s on the right side due to higher CD8). Runs downstream stats based on ratio.

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(ggsignif)
library(ggpmisc)

##..Data importation, clean-up and reshaping..## 

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

data<-read.csv("immune_cell_freqs.csv")

#reshape matrix so that each row is a point, each column is a frequency for that cell type
freq_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "frequency") #frequency data
count_data <- reshape2::dcast(data, SampleID ~ cell_type, value.var = "count") #count data

#create T cell count only matrix 
tcell_data<-count_data[,c(1,10,12)]

#add sample information to all matrices
Samp_num<-c(30,30,11,17,17,23,2,20,20,4,18,18,34,34,31,3,29,29,12,12,21,21,31,1,1,67,68,69,70,71,72,73,74,75,76,2,11,23,3,4)
tcell_data$samp_num<-Samp_num
freq_data$samp_num<-Samp_num
count_data$samp_num<-Samp_num

#annotate origin of sample based on sample number
AHRI<-c(1,2,3,4,11,34)
Stanford<-c(17,21,23,30,18,20,31,29,12)
Sarcoid<-c(67,68,69,70,71,72,73,74,75,76)

tcell_data<-tcell_data %>% mutate(Origin=case_when(tcell_data$samp_num %in% AHRI ~ "AHRI",
                                                   tcell_data$samp_num %in% Stanford ~ "Stanford",
                                                   tcell_data$samp_num %in% Sarcoid ~ "Sarcoid"))
freq_data<-freq_data %>% mutate(Origin=case_when(freq_data$samp_num %in% AHRI ~ "AHRI",
                                                  freq_data$samp_num %in% Stanford ~ "Stanford",
                                                 freq_data$samp_num %in% Sarcoid ~ "Sarcoid"))
count_data<-count_data %>% mutate(Origin=case_when(count_data$samp_num %in% AHRI ~ "AHRI",
                                                   count_data$samp_num %in% Stanford ~ "Stanford",
                                                   count_data$samp_num %in% Sarcoid ~ "Sarcoid"))

##..Calculate log2(CD4/CD8)..##

log2fc<-log2(tcell_data$CD4_T/tcell_data$CD8_T)
tcell_data$FC<-log2fc

##..Plot ratio and color by sample of origin or sample site..##

# TB only
MAC_IDs<-c(57,58,19,87)
tcell_data_TB<-droplevels(tcell_data[!tcell_data$Origin=="Sarcoid",])
tcell_data_TB<-droplevels(tcell_data_TB[!tcell_data_TB$SampleID %in% MAC_IDs,])

# TB and Sarcoid
tcell_data_noMAC<-droplevels(tcell_data[!tcell_data$SampleID %in% MAC_IDs,])
tcell_data_noMAC<-tcell_data_noMAC %>% mutate(Tissue=case_when(tcell_data_noMAC$samp_num %in% AHRI ~ "TB",
                                                               tcell_data_noMAC$samp_num %in% Stanford ~ "TB",
                                                               tcell_data_noMAC$samp_num %in% Sarcoid ~ "Sarcoid"))

# Append cluster to the TB data for visualization

cluster_data<-as.data.frame(tcell_data_TB$SampleID)
colnames(cluster_data)<-'SampleID'
cluster1<-c(33,34,47,48,6,7,61,26,27)
cluster2<-c(15,21,40,84)
cluster3<-c(14,42,88,13,89,85,64,65,55,28,36,35,54)
cluster_data$cluster<-NA
cluster_data[cluster_data$SampleID %in% cluster1,]$cluster<-1
cluster_data[cluster_data$SampleID %in% cluster2,]$cluster<-2
cluster_data[cluster_data$SampleID %in% cluster3,]$cluster<-3
tcell_data_TB$cluster<-cluster_data$cluster

# Assign color palettes

color_disease<-c("#00A59C","#9CD9D3") #TB and sarcoid
color_specimen<-c('#67AF85','#5982AF') #resection and biopsy
color_cluster<-c('#2E3192','#006838','#9e1f63') #cluster

# T cell ratio in all TB FOVs in descending order
tcell_ratio<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -FC), y=FC)) +
  geom_bar(stat="Identity", aes(fill=as.factor(cluster))) +
  scale_fill_manual(values=color_cluster) +
  geom_hline(yintercept=c(0,0), linetype="dashed") +
  theme_bw() + 
  theme(legend.position = 'none') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="CD4:CD8 Ratio (log2 [CD4 T cells / CD8 T cells])")
tcell_ratio

# CD4:CD8 Ratio in Sarcoid v TB
ratio_sarcvTB<-ggplot(data=tcell_data_noMAC, aes(x=Tissue, y=FC, fill=Tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 2, position = position_jitterdodge()) +
  theme_bw() + 
  scale_fill_manual(values=color_disease) +
  stat_compare_means(data=tcell_data_noMAC,label = "p.signif", method= "wilcox.test") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(legend.position = 'none') +
  labs(y="FC Sarcoid v TB")
ratio_sarcvTB

##..Plot and compare the absolute CD4 and CD8 counts as well as the frequencies..##ses s

# Total CD4 counts per sample

CD4_count_persample<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -CD4_T), y=CD4_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD4 T cells") 
CD4_count_persample


# Total CD8 counts per sample

CD8_count_persample<-ggplot(data=tcell_data_TB, aes(x=reorder(SampleID, -CD8_T), y=CD8_T, fill=Origin)) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_specimen) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total Number of CD8 T cells")
CD8_count_persample


## Frequency ##
FC_order<-c(15,27,55,6,33,34,26,61,40,7,47,88,42,84,21,28,54,14,85,13,36,35,89,65,48,64)
freq_data_TB<-droplevels(freq_data[freq_data$samp_num %in% c(Stanford,AHRI),])
freq_data_TB<-droplevels(freq_data_TB[!freq_data_TB$SampleID %in% MAC_IDs,])
freq_data_TB$cluster<-tcell_data_TB$cluster
freq_data_TB$SampleID<-factor(freq_data_TB$SampleID, levels = FC_order)
freq_data_TB<-freq_data_TB[order(freq_data_TB$SampleID),]

# CD4 frequency per sample

CD4_freq_persample<-ggplot(data=freq_data_TB, aes(x=as.factor(SampleID), y=CD4_T, fill=as.factor(cluster))) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_cluster) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)") +
  theme(legend.position = 'none')
CD4_freq_persample

# CD8 frequency per sample

CD8_freq_persample<-ggplot(data=freq_data_TB, aes(x=as.factor(SampleID), y=CD8_T, fill=as.factor(cluster))) +
  geom_bar(stat="Identity") +
  theme_bw() + 
  scale_fill_manual(values=color_cluster) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)") +
  theme(legend.position = 'none')
CD8_freq_persample


# Cluster comparisons

data_summary <- freq_data_TB %>%
  group_by(cluster) %>%
  summarize(combo_median = median(CD4_T))

my_comparisons = list(c(1,2),
                   c(2,3),
                   c(3,1))

CD4_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(cluster), y=CD4_T, color=as.factor(cluster))) +
  geom_point(data=freq_data_TB, aes(x=as.factor(cluster), y=CD4_T, color=as.factor(cluster)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_cluster) +
  geom_point(aes(y = combo_median), color = "black", size = 2, data =data_summary) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD4 T cells (of total immune cells)")
CD4_freq_bulk


data_summary <- freq_data_TB %>%
  group_by(cluster) %>%
  summarize(combo_median = median(CD8_T))


CD8_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(cluster), y=CD8_T, color=as.factor(cluster))) +
  geom_point(data=freq_data_TB, aes(x=as.factor(cluster), y=CD8_T, color=as.factor(cluster)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_cluster) +
  geom_point(aes(y = combo_median), color = "black", size = 2, data =data_summary) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD8 T cells (of total immune cells)")
CD8_freq_bulk


# CD14 monos

data_summary <- freq_data_TB %>%
  group_by(cluster) %>%
  summarize(combo_median = median(CD14_Mono))


CD14_Mono_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(cluster), y=CD14_Mono, color=as.factor(cluster))) +
  geom_point(data=freq_data_TB, aes(x=as.factor(cluster), y=CD14_Mono, color=as.factor(cluster)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_cluster) +
  geom_point(aes(y = combo_median), color = "black", size = 2, data =data_summary) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of CD14+ Monos (of total immune cells)")
CD14_Mono_freq_bulk


data_summary <- freq_data_TB %>%
  group_by(cluster) %>%
  summarize(combo_median = median(`CD11b/c_CD206_Mac/Mono`))

MDSC_freq_bulk<-ggplot(data=freq_data_TB, aes(x=as.factor(cluster), y=`CD11b/c_CD206_Mac/Mono`, color=as.factor(cluster))) +
  geom_point(data=freq_data_TB, aes(x=as.factor(cluster), y=`CD11b/c_CD206_Mac/Mono`, color=as.factor(cluster)), 
             position = position_jitter(width = 0.1, height = 0.0), size = 6) + 
  scale_color_manual(values=color_cluster) +
  geom_point(aes(y = combo_median), color = "black", size = 2, data =data_summary) +
  geom_errorbar(aes(y = combo_median, ymin = combo_median, ymax = combo_median),
                color = "black", width = 0.2, data = data_summary) +
  theme_bw() + 
  theme(legend.position = 'none') +
  stat_compare_means(label = "p.format", method= "wilcox.test", comparisons = my_comparisons) +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Frequency of MDSC (of total immune cells)")
MDSC_freq_bulk



# AHRI vs Stanford CD4 count

CD4_count_bulk<-ggplot(data=tcell_data_TB, aes(x=Origin, y=CD4_T, color=Origin)) +
  geom_boxplot(lwd=1) +
  theme_bw() + 
  stat_compare_means(label = "p.format",method= "wilcox.test", comparisons = list(c("AHRI","Stanford"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total CD4 T cells")
CD4_count_bulk

# AHRI vs Stanford CD4 freq

CD4_freq_bulk<-ggplot(data=freq_data_TB, aes(x=Origin, y=CD4_T, color=Origin)) +
  geom_boxplot(lwd=1) +
  theme_bw() + 
  stat_compare_means(label = "p.format",method= "wilcox.test", comparisons = list(c("AHRI","Stanford"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Freq CD4 T cells")
CD4_freq_bulk

# AHRI vs Stanford CD8 count

CD8_count_bulk<-ggplot(data=tcell_data_TB, aes(x=Origin, y=CD8_T, color=Origin)) +
  geom_boxplot(lwd=1) +
  theme_bw() + 
  stat_compare_means(label = "p.format",method= "wilcox.test", comparisons = list(c("AHRI","Stanford"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Total CD8 T cells")
CD8_count_bulk

# AHRI vs Stanford CD8 freq

CD8_freq_bulk<-ggplot(data=freq_data_TB, aes(x=Origin, y=CD8_T, color=Origin)) +
  geom_boxplot(lwd=1) +
  theme_bw() + 
  stat_compare_means(label = "p.format",method= "wilcox.test", comparisons = list(c("AHRI","Stanford"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x=element_blank(), 
        axis.title.x=element_blank(), axis.ticks.x=element_blank()) + 
  labs(y="Freq CD8 T cells")
CD8_freq_bulk

# CD4/CD8 Ratio versus the CD4 or CD8 Count (as a scatter plot with linear regression and pearson corr)

fitCD4 <- lm(CD4_T ~ FC, data = tcell_data)
corr_CD4<-cor.test(tcell_data$CD4_T,tcell_data$FC,method="pearson")
corrCoeff_CD4<-corr_CD4$estimate

ggplot(tcell_data,aes(x=CD4_T,y=FC)) + geom_point(aes(color=Origin),size=3) + labs(x="CD4 T Cell Count") + 
  geom_smooth(method='lm',formula=y~x,alpha=0.25,colour="black") +
  labs(y="CD4/CD8 Ratio") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title = paste("Adj R2 = ",signif(summary(fitCD4)$adj.r.squared, 2),
                     " P =",signif(summary(fitCD4)$coef[2,4], 2), 
                     " Pearson R =",as.numeric(round(corrCoeff_CD4,digits=2))))


fitCD8 <- lm(FC ~ CD8_T, data = tcell_data)
corr_CD8<-cor.test(tcell_data$CD8_T,tcell_data$FC,method="pearson")
corrCoeff_CD8<-corr_CD8$estimate

ggplot(tcell_data,aes(x=CD8_T,y=FC)) + geom_point(aes(color=Origin),size=3) + labs(x="CD8 T Cell Count") + 
  geom_smooth(method='lm',formula=y~x,alpha=0.25,colour="black") +
  labs(y="CD4/CD8 Ratio") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title = paste("Adj R2 = ",signif(summary(fitCD8)$adj.r.squared, 2),
                     " P =",signif(summary(fitCD8)$coef[2,4], 2), 
                     " Pearson R =",as.numeric(round(corrCoeff_CD8,digits=2))))


# CD4 versus Treg Counts
fitTreg <- lm(Treg ~ CD4_T, data = count_data)
corr_Treg<-cor.test(count_data$CD4_T,count_data$Treg,method="pearson")
corrCoeff_Treg<-corr_Treg$estimate

ggplot(count_data,aes(x=Treg,y=CD4_T)) + geom_point(aes(color=Origin),size=3) + labs(x="Treg Count") + 
  geom_smooth(method='lm',formula=y~x,alpha=0.25,colour="black") +
  labs(y="CD4 T Cell Count") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title = paste("Adj R2 = ",signif(summary(fitTreg)$adj.r.squared, 2),
                     " P =",signif(summary(fitTreg)$coef[2,4], 2), 
                     " Pearson R =",as.numeric(round(corrCoeff_Treg,digits=2))))

# CD4 versus Treg Freq
fitTreg2 <- lm(Treg ~ CD4_T, data = freq_data)
corr_Treg2<-cor.test(freq_data$CD4_T,freq_data$Treg,method="pearson")
corrCoeff_Treg2<-corr_Treg2$estimate

ggplot(freq_data,aes(x=Treg,y=CD4_T)) + geom_point(aes(color=Origin),size=3) + labs(x="Treg Frequency") + 
  geom_smooth(method='lm',formula=y~x,alpha=0.25,colour="black") +
  labs(y="CD4 T Cell Frequency") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  labs(title = paste("Adj R2 = ",signif(summary(fitTreg2)$adj.r.squared, 2),
                     " P =",signif(summary(fitTreg2)$coef[2,4], 2), 
                     " Pearson R =",as.numeric(round(corrCoeff_Treg2,digits=2))))


# Correlation between CD4 and CD8 T cells

##..Get rank sum p-value between AHRI and Stanford..##

tcell_data_noMAC<-tcell_data_noMAC %>% mutate(Tissue=case_when(tcell_data_noMAC$samp_num %in% AHRI ~ "TB",
                                                               tcell_data_noMAC$samp_num %in% Stanford ~ "TB",
                                                               tcell_data_noMAC$samp_num %in% Sarcoid ~ "Sarcoid"))

compare_means(FC ~ Origin, data = tcell_data_noMAC[tcell_data_noMAC$Tissue == 'TB',], method = 'wilcox.test')

##..Save dataframe SampleID, Samp_num, CD4 total, CD8 total, CD4/CD8 ratio and save..##

write.csv(tcell_data, file="CD4CD8ratio.csv",row.names = FALSE)
