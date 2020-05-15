# TBMIBIcellCountChiSquare.R
# Author: Erin McCaffrey 
# Date created: 200419
# Overview: This script reads in the fcount data for all TB granulomas and runs a chi-square
# between all pairs of cell types to assess the significance of co-occurance.

library(reshape2)
library(dplyr)
library(tidyverse)
library(broom)
library(gplots)
library(RColorBrewer)

# import data
setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")

#cell frequencies
imm_freqs<-read.csv("immune_cell_freqs.csv")
lin_freqs<-read.csv("lineage_freqs.csv")

#reshape to a matrix where rows are samples and columns are counts of cell types
immfreq_data <- reshape2::dcast(imm_freqs, SampleID ~ cell_type, value.var = "count")
linfreq_data <- reshape2::dcast(lin_freqs, SampleID ~ lineage, value.var = "count")

#merge into one matrix with epithelial, endothelial, and fibroblast data
combo<-merge(immfreq_data,linfreq_data[,-5],by="SampleID")

#remove MAC samples
MAC_IDs<-c(57,58,19,87)
combo<-combo[!combo$SampleID %in% MAC_IDs,]

#drop sarcoid
sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
combo <- combo[!combo$SampleID %in% sarcoid,]

#drop SampleID
combo<-combo[,-1]

#produce version based on presence/absence
t<-5 #number of cells required for positivity
combo_binary<-combo
combo_binary<-ifelse(combo_binary>=5,1,0)

# run chi square over all column pairs

df<-combo
list.a<-colnames(df)
list.b<-colnames(df)
combos<-expand.grid(list.a, list.b)
colnames(combos)<-c('X1','X2')
combos<-subset(combos, X1 != X2)
  
chisquare_allpairs<-combos %>%
  mutate(d = map2(X1, X2, ~tidy(chisq.test(df[,.x], df[,.y])))) %>%
  unnest()

# sort out unique combinations

chisquare_uniquepairs <- chisquare_allpairs %>%
  mutate(Var = map2_chr(X1, X2, ~toString(sort(c(.x, .y))))) %>%
  distinct(Var, .keep_all = TRUE) %>%
  select(-Var)

# get adjusted p-values
chisquare_uniquepairs$adjp<-p.adjust(chisquare_uniquepairs$p.value, method='fdr')

# cast to a grid to plot as a matrix
chisq_grid <- reshape2::acast(chisquare_uniquepairs, X1 ~ X2, value.var = "adjp")


# plot with p-values and then binary sig or not
row_order<-row.names(chisq_grid)
colorkey<-read.csv('colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% row_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = row_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color_row<-as.vector(colorkey_imm$code)

col_order<-colnames(chisq_grid)
colorkey<-read.csv('colorkey_R.csv')
colorkey_imm<-droplevels(colorkey[colorkey$imm_order %in% col_order,])
colorkey_imm$imm_order<-factor(colorkey_imm$imm_order, levels = col_order)
colorkey_imm<-colorkey_imm[order(colorkey_imm$imm_order),]
color_col<-as.vector(colorkey_imm$code)


heatmap.2(chisq_grid,
          col = colorRampPalette((brewer.pal(9,"Blues")))(100),,
          trace = "none",
          Colv = F, Rowv = F,
          dendrogram = "none",
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=0:ncol(chisq_grid),
          rowsep=0:nrow(chisq_grid),
          cellnote = chisq_grid,
          notecex=0.3,
          notecol="black",
          na.color=par("bg"),
          cexRow=0.4,
          cexCol=0.4,
          RowSideColors = color_row,
          ColSideColors = color_col,
          symm = TRUE)

#convert to binary matrix

chisq_grid_binary<-chisq_grid
chisq_grid_binary[chisq_grid_binary==0.000]<-1
chisq_grid_binary<-ifelse(chisq_grid_binary<0.05,1,0)

heatmap.2(chisq_grid_binary,
          col = colorRampPalette((brewer.pal(9,"Blues")))(100),,
          trace = "none",
          Colv = F, Rowv = F,
          dendrogram = "none",
          density.info = 'none',
          key.title = '',
          sepwidth=c(0.01,0.01),
          sepcolor="black",
          colsep=0:ncol(chisq_grid_binary),
          rowsep=0:nrow(chisq_grid_binary),
          cexRow=0.4,
          cexCol=0.4,
          RowSideColors = color_row,
          ColSideColors = color_col,
          symm = TRUE)


