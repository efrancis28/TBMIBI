# TBMIBIrunUMAP.R
# Author: Erin McCaffrey
# Date created: 190116
# Overview: This script reads in the csv for cell-size normalized, asinh transformed, and percentile normalized data.
# Next it runs UMAP on chosen data and visualizes results colored by cell type, sample ID, or tissue type.

library("umap")
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library("ggrastr")

##..Import data..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data<-read.csv("granA_cellpheno_CS-asinh-norm_revised.csv")
# mycobacterial samples only
data_myco<-droplevels(data[data$Tissue %in% c('gran_lung','gran_pleura','gran_endo','gran_LN'),])


##..Choose clustering channels for different implementations..##

clusterChannels = c('CD14','CD16','CD68','CD163','CD11b','CD11c','CD209','CD206','CD45')

#myeloid cells#

# get just myeloid cells #
myeloid<-c("CD14_Mono","CD11b/c_CD206_Mac/Mono","CD11c_DC/Mono","CD68_Mac","CD16_CD14_Mono",
           "CD206_Mac","CD163_Mac","CD209_DC","giant_cell")
data_myeloid<-droplevels(data_myco[data_myco$cell_type %in% myeloid,])

set.seed(321)
gran_myeloid.umap<-umap(data_myeloid[,clusterChannels])
umap.output2<-gran_myeloid.umap$layout
data_myeloid$umap1<-umap.output2[,1]
data_myeloid$umap2<-umap.output2[,2]

##..Visualize..##

#write.csv(data_myeloid, file="granA_myeloid_umap.csv",row.names = FALSE)
data_myeloid<-read.csv("granA_myeloid_umap.csv")

#Color by sample ID, tissue, or cell ID
data_myeloid$cell_type<-factor(data_myeloid$cell_type, levels=myeloid)
data_myeloid<-data_myeloid[order(data_myeloid$cell_type),]

##..Plot..##
colorkey<-read.csv("colorkey_R.csv")
myeloid_color<-droplevels(colorkey[colorkey$imm_order %in% myeloid,])
myeloid_color$imm_order<-factor(myeloid_color$imm_order, levels = myeloid)
color<-as.vector(myeloid_color$code)


ggplot(data_myeloid, aes(umap1, umap2, color = as.factor(cell_type))) +
  geom_point(size=2.75)  + 
  theme_bw() + 
  scale_color_manual(values = color) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  labs(x="UMAP 1") +
  labs(y="UMAP 2") +
  labs(color = "Immune Phenotype") +
  ggtitle("UMAP of Mycobacterial Myeloid Cells")


#overlay marker expression 

theme <- theme(strip.background = element_blank(),
               panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank())

data_myeloid.melt <- reshape2::melt(data_myeloid, id.vars = c('umap1', 'umap2'), measure.vars = c(clusterChannels2,
                                                                                                  'IDO','PD.L1','HLA.DR.DQ.DP','Ki67','H3K9Ac','CD36' ,'iNOS','Collagen.1'))

ggplot(data_myeloid.melt, aes(x = umap1, y = umap2, color = value)) +
  geom_point_rast(size=3) +
  scale_color_viridis(option = 'magma', limits=c(0, 1)) +
  facet_wrap(~ variable, ncol = 4) +
  theme +
  theme(legend.position = 'right')

