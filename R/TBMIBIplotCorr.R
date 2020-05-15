# TBMIBIplotCorr.R
# Author: Erin McCaffrey 
# Date created: 190319
# Overview: This script reads in the csv for cell-size normalized data and then plots two markers, 
# determines their Pearson correlation and fits with a linear regression. 

library(ggplot2)

#read in data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")

#keep mycobacterial samples and myeloid cells only

data_myco<-droplevels(data[data$Tissue  %in% c('gran_lung','gran_pleura','gran_endo','gran_LN'),])
myeloid<-c("CD14_Mono","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","CD16_CD14_Mono",
           "CD206_Mac","CD163_Mac","CD209_DC","giant_cell")
data_myeloid<-droplevels(data_myco[data_myco$cell_type %in% myeloid, ])

##..Run pearson correlation and get coefficient..##

corr<-cor.test(data_myco$IDO,data_myco$PD.L1,method="pearson")
corrCoeff<-corr$estimate

corr_myeloid<-cor.test(data_myeloid$IDO,data_myeloid$PD.L1,method="pearson")
corrCoeff_myeloid<-corr_myeloid$estimate

##..Break down correlation by cohort and tissue..##

resection<-c(64,65,21,84,42,88,28,89,85,13,35,36)
corr_resection<-cor.test(data_myco[data_myco$SampleID %in% resection,]$IDO,
                         data_myco[data_myco$SampleID %in% resection,]$PD.L1,method="pearson")
corrCoeff_resection<-corr_resection$estimate

diagnostic_pulm<-c(6,7,14,15)
corr_diagnostic_pulm<-cor.test(data_myco[data_myco$SampleID %in% diagnostic_pulm,]$IDO,
                         data_myco[data_myco$SampleID %in% diagnostic_pulm,]$PD.L1,method="pearson")
corrCoeff_diagnostic_pulm<-corr_diagnostic_pulm$estimate

pulm<-c(resection, diagnostic_pulm)
corr_pulm<-cor.test(data_myco[data_myco$SampleID %in% pulm,]$IDO,
                         data_myco[data_myco$SampleID %in% pulm,]$PD.L1,method="pearson")
corrCoeff_pulm<-corr_pulm$estimate

diagnostic_expulm<-c(33,34,26,27,40,61,47,48,54,55)
corr_diagnostic_expulm<-cor.test(data_myco[data_myco$SampleID %in% diagnostic_expulm,]$IDO,
                         data_myco[data_myco$SampleID %in% diagnostic_expulm,]$PD.L1,method="pearson")
corrCoeff_diagnostic_expulm<-corr_diagnostic_expulm$estimate

diagnostic<-c(diagnostic_pulm, diagnostic_expulm)
corr_diagnostic<-cor.test(data_myco[data_myco$SampleID %in% diagnostic,]$IDO,
                                 data_myco[data_myco$SampleID %in% diagnostic,]$PD.L1,method="pearson")
corrCoeff_diagnostic<-corr_diagnostic$estimate

CD4skewed<-c(6,7,15,26,27,33,34,40,55,6)
corr_CD4<-cor.test(data_myco[data_myco$SampleID %in% CD4skewed,]$IDO,
                          data_myco[data_myco$SampleID %in% CD4skewed,]$PD.L1,method="pearson")
corrCoeff_CD4<-corr_CD4$estimate

CD8skewed<-c(13,14,35,36,48,54,64,65,85,89)
corr_CD8<-cor.test(data_myco[data_myco$SampleID %in% CD8skewed,]$IDO,
                          data_myco[data_myco$SampleID %in% CD8skewed,]$PD.L1,method="pearson")
corrCoeff_CD8<-corr_CD8$estimate


##..Run pearson correlation on CD4 and CD8 skewed group separately and get coefficient..##

skew_data<-read.csv("CD4CD8ratio.csv")
TB_samples <- unique(data_myco$SampleID)
skew_data<-droplevels(skew_data[skew_data$SampleID %in% TB_samples, ])
CD4skewed<-skew_data[which(skew_data$FC>0.76),]$SampleID
CD8skewed<-skew_data[which(skew_data$FC<(-0.40)),]$SampleID

skew_data_exprs <- droplevels(data_myeloid[data_myeloid$SampleID %in% c(CD4skewed, CD8skewed), ])
skew_data_exprs$TcellRatio <- "CD4_skewed"
skew_data_exprs[skew_data_exprs$SampleID %in% CD8skewed, ]$TcellRatio <- "CD8_skewed"

corr_myeloid_CD8 <-cor.test(data_myeloid[data_myeloid$SampleID %in% CD8skewed,]$IDO,data_myeloid[data_myeloid$SampleID %in% CD8skewed,]$PD.L1,method="pearson")
corrCoeff_myeloid_CD8<-corr_myeloid_CD8$estimate

corr_myeloid_CD4 <-cor.test(data_myeloid[data_myeloid$SampleID %in% CD4skewed,]$IDO,data_myeloid[data_myeloid$SampleID %in% CD4skewed,]$PD.L1,method="pearson")
corrCoeff_myeloid_CD4<-corr_myeloid_CD4$estimate

##..Plot IDO and PD.L1 with correlation and linear regression..##

corrCoeff_txt<-as.numeric(round(corrCoeff,digits=2))
corrString=toString(corrCoeff_txt)


ggplot(data_myeloid,aes(x=IDO,y=PD.L1)) +
  geom_point() + 
  labs(x="IDO Expression") + 
  labs(y="PD-L1 Expression") + theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 


