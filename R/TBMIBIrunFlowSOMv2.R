# TBMIBIrunFlowSOM.R
# Author: Erin McCaffrey (with many sections adapted from Felix's R Demo)
# Date created: 190128
# Overview: Script imports concatenated data set of asinh transformed sc data. Runs FlowSOM on the transformed and
# normalized data.


                ###########################################
                ##..Install packages and open libraries..##
                ###########################################

#install.packages("BiocManager")
#BiocManager::install("BiocUpgrade")  # upgrading Bioconductor
#BiocManager::install("flowCore")     # for handling of fcs files in R
#BiocManager::install("ggplot2")      # for advanced data plotting
#BiocManager::install("gplots")       # for heatmaps
#BiocManager::install("RColorBrewer") # additional color schemes
#BiocManager::install("reshape2")     # reorganizing data
#BiocManager::install("FlowSOM")      # FlowSOM clustering
#BiocManager::install("plyr")         # for data handling
#install.packages("viridis")          # for badass color palettes

##..Open necessary libraries..##

library(flowCore)           
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(FlowSOM)
library(plyr)
library(dplyr)
library(viridis)

                  ################################################################################
                  ##..Data Preparation of A granulomas, CS-normalized, scaled, and asinh trans..##
                  ################################################################################

##..Read in the concatenated expression matrix..##

setwd("/Volumes/GoogleDrive/My Drive/Angelo Lab/MIBIProjects/Human_ATB_paper-cohort/all-fields/denoised/dataPerCell_3px")
data_CSscaled<-read.csv("allsamples_dataCS3px_annotated.csv") #cell size normalized and untransformed

##..Include only A granulomas..##

granA_lung<-c(64,65,21,84,42,88,28,89,85,13,35,36,14,15,57,58,19,87,6,7)
granA_pleura<-c(33,34,26,27,40,61)
granA_endo<-c(47,48)
granA_LN<-c(54,55)
gran_sarcoid<-c(67,68,69,70,71,72,73,74,75,76)
all<-c(granA_lung,granA_pleura,granA_LN,granA_endo,gran_sarcoid) #all A granulomas

data_gran<-data_CSscaled[data_CSscaled$SampleID %in% all,] #exprs data only A granulomas
#write.csv(data_gran, file="gran_dataCS_scale.csv",row.names = FALSE)

# linear transform the data by 100 (only the expression data and not sample ID, cell label, cell size, or tissue)

data_gran[,4:50]<- data_gran[,4:50]*100

# arcsinh transform data (only the expression data and not sample ID, cell label, cell size, or tissue)

asinh_scale <- 5
data_trans<-data_gran
data_trans[,4:50]<- asinh(data_gran[,4:50]/ asinh_scale)

##..Percent normalization and scale data 0-99th percentile..##

v <- 1:1000
v
quantile_value <- 0.999
quantile(v, quantile_value)

# calculating percentile for 1 vector
percentile.vector <- apply(data_trans[,4:50], 2, function(x) quantile(x, quantile_value, names = F))
percentile.vector
data_gran_norm<-data_trans
data_gran_norm[,4:50] <- data.frame(t(t(data_trans[,4:50]) / as.numeric(percentile.vector)))
data_gran_norm<-data_gran_norm[,-c(7,8,17,34)] #remove channels with no expression (Ca, Fe, arg1, CD56)



                #######################################################################
                ###..FlowSOM Cluster round 1: Immune, Endo, Epithelial, Fibroblast..###
                #######################################################################

ff_new <- flowFrame(exprs = data.matrix(data_gran_norm[,-c(47,48)]), desc = list(FIL = 1)) #exclude tissue type column

clusterChannels_1<-c('CD45','SMA','Vimentin','CD31','Keratin.pan','E.cadherin','MastChyTry','MPO','CD20','CD3','CD14','HLA.DR.DQ.DP')

##..Run FlowSOM random seed for reproducibility..##

set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clusterChannels_1)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
fs_clusters <- out_fSOM$map$mapping[,1]

out_fSOM <- UpdateNodeSize(out_fSOM, reset = TRUE)
FlowSOM::PlotStars(out_fSOM, view = "grid", markers = clusterChannels_1)
out_fSOM <- UpdateNodeSize(out_fSOM)
FlowSOM::PlotStars(out_fSOM, view = "MST", markers = clusterChannels_1)
FlowSOM::PlotStars(out_fSOM, view="tSNE",markers = clusterChannels_1)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
data_fs_clusters <- cbind(data_gran_norm, fs_clusters)

# go through all clusters and calculate mean for every channel
hm_allclusters <- matrix(, nrow = length(unique(fs_clusters)), ncol = length(clusterChannels_1))
for(i in 1:length(unique(fs_clusters))) {
  temp_mat <- data_fs_clusters[data_fs_clusters[,"fs_clusters"] == i, clusterChannels_1]
  hm_allclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}
hm_allclusters[is.na(hm_allclusters)] <- 0

# add names to rows and cols
rownames(hm_allclusters) <- paste("cluster", 1:length(unique(fs_clusters)), sep = "")
colnames(hm_allclusters) <- clusterChannels_1  
hm_allclusters

# plot heatmap of all clusters
#tiff("plots/draft_figs/190605_FlowSOM-lineage-hmap.tiff", units="in", width=15, height=11, res=300)
#setEPS()
#postscript("plots/draft_figs/190605_FlowSOM-lineage-hmap.eps",width=15, height=11)
heatmap.2(hm_allclusters, 
           scale = "none",
           Colv = T, Rowv = T,
           hclustfun = hclust,
           dendrogram = c("both","row","column","none"),
           trace = "none",
           #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
           #col = colorRampPalette(brewer.pal(9,"Blues"))(100),
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.3, cexCol = 0.4, 
          breaks=seq(0, 1, length.out=257))
#dev.off()

# Pull out sparse pops just in case they are absorbed in clustering

neut_row_idx<-which(data_fs_clusters$fs_clusters == 45)
neut_rows<-rownames(data_fs_clusters[neut_row_idx,])

mast_row_idx<-which(data_fs_clusters$fs_clusters == 56)
mast_rows<-rownames(data_fs_clusters[mast_row_idx,])


##..Meta-cluster..##

# try the suggested automatic metaclustering method for a hint for k
auto_meta <- MetaClustering(out_fSOM$map$codes, method = "metaClustering_consensus", max = 20)
max(auto_meta)

# do a manual metaclustering
chosen_k=20
set.seed(123)
out_meta <- FlowSOM::metaClustering_consensus(out_fSOM$map$codes, k = chosen_k)
meta_results <- out_meta[fs_clusters]

##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters <- cbind(data_gran_norm, meta_results)

# check the amount of cells in each cluster
table(data_meta_clusters$meta_results)

# go through all clusters and calculate mean for every channel
hm_metaclusters <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_1))
for(i in 1:chosen_k) {
   temp_mat <- data_meta_clusters[data_meta_clusters[,"meta_results"] == i, clusterChannels_1]
   hm_metaclusters[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
 }

# rename
rownames(hm_metaclusters) <- paste("cluster", 1:chosen_k, sep = "")
colnames(hm_metaclusters) <- clusterChannels_1
hm_metaclusters

# make a metacluster heatmap

color<-c(rep('#3366FF',3),'#99FF99',rep('#3366FF',8),'#CC3333','#3366FF','#FF9966',
         '#99FF99','#FF9966',rep('#CC3333',3))

heatmap.2(hm_metaclusters,
        scale = "none",
        Colv = T, Rowv = T,
        hclustfun = hclust,
        dendrogram = c("both"),
        trace = "none",
        col = viridis(256),
        sepwidth=c(0.01,0.01),
        sepcolor="grey",
        colsep=0:ncol(hm_allclusters),
        rowsep=0:nrow(hm_allclusters),
        RowSideColors=color,
        density.info = 'none',
        key.title = '',
        cexRow = 0.7, cexCol = 0.7, 
        breaks=seq(0, 1, length.out=257))


##..Plot the percent of all cells in each major category..##

imm_clust<-c(1,2,3,5,6,7,8,9,10,11,12,14)
fibro_clust<-c(4,16)
endo_clust<-c(13,18,19,20)
epi_clust<-c(15,17)
cell_type<-meta_results

cell_type<-replace(cell_type,cell_type %in% imm_clust,"immune")
cell_type<-replace(cell_type,cell_type %in% epi_clust,"epithelial")
cell_type<-replace(cell_type,cell_type %in% endo_clust,"endothelial")
cell_type<-replace(cell_type,cell_type %in% fibro_clust,"fibroblast")
freq<-data.frame(table(cell_type))

cell_freq<-freq
cell_freq$percent<-as.numeric(format((cell_freq$Freq/sum(cell_freq$Freq)*100), digits=3))

ggplot(data=cell_freq, aes(x=reorder(cell_type, -percent), y=percent, fill=cell_type)) +
  geom_text(aes(label=percent), vjust=-0.3, size=3.5)+
  geom_bar(stat="identity") +
  scale_fill_manual(values =c('#CC3333','#FF9966','#99FF99','#3366FF')) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Cell Type") +
  labs(y="Percent of Total") +
  ggtitle("Frequency of Cell Types") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y=c(0,100))

##..Add cell type to each event..##
data_gran$lineage<-cell_type
data_gran_norm$lineage<-data_gran$lineage

  
                ###################################################
                ##..FlowSOM Cluster round 2: Major Immune cells..##
                ###################################################

library(dplyr)
immune<-c("immune") #can append if wanting to include non-immune clusters that look contaminated
data_immune<-data_gran_norm[data_gran_norm$lineage %in% immune,]
ff_immune <- flowFrame(exprs = data.matrix(data_immune[,-c(47,48,49)]), desc = list(FIL = 1))

clusterChannels_2=c("CD3","CD20","CD4","CD8","Foxp3","CD14","CD16","CD11c","MastChyTry","MPO","CD206")

##..Run FlowSOM random seed for reproducibility..##

set.seed(123)
out_fSOM_imm <- FlowSOM::ReadInput(ff_immune, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_imm <- FlowSOM::BuildSOM(out_fSOM_imm, colsToUse = clusterChannels_2)
out_fSOM_imm <- FlowSOM::BuildMST(out_fSOM_imm)
labels_imm <- out_fSOM_imm$map$mapping[,1]

out_fSOM_imm <- UpdateNodeSize(out_fSOM_imm, reset = TRUE)
FlowSOM::PlotStars(out_fSOM_imm, view = "grid", markers = clusterChannels_2)
out_fSOM_imm <- UpdateNodeSize(out_fSOM_imm)
FlowSOM::PlotStars(out_fSOM_imm, view = "MST", markers = clusterChannels_2)
FlowSOM::PlotStars(out_fSOM_imm, view="tSNE",markers = clusterChannels_2)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_imm<-out_fSOM_imm[["map"]][["mapping"]]
fs_clusters_imm<-fs_clusters_imm[,1]
data_fs_clusters_imm <- cbind(data_immune, fs_clusters_imm)

# go through all clusters and calculate mean for every channel
hm_allclusters_imm <- matrix(, nrow = length(unique(fs_clusters_imm)), ncol = length(clusterChannels_2))
for(i in 1:length(unique(fs_clusters_imm))) {
  temp_mat <- data_fs_clusters_imm[data_fs_clusters_imm[,"fs_clusters_imm"] == i, clusterChannels_2]
  hm_allclusters_imm[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_imm) <- paste("cluster", 1:length(unique(fs_clusters_imm)), sep = "")
colnames(hm_allclusters_imm) <- clusterChannels_2  
hm_allclusters_imm

# plot heatmap of all clusters

imm_hmap_color<-sort(unique(fs_clusters_imm))
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% myeloid_clust,"#808080")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% CD4_clust,"#00B3F2")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% CD8_clust,"#45B500")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% B_clust,"#00C087")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% treg_clust,"#9C8DFF")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% mast_clust,"#FF61C7")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% neut_clust,"#EBBC00")
imm_hmap_color<-replace(imm_hmap_color,imm_hmap_color %in% other_clust,"#B2A100")


heatmap.2(hm_allclusters_imm, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both"),
          trace = "none",
          col = viridis(256),
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_allclusters_imm),
          rowsep=0:nrow(hm_allclusters_imm),
          RowSideColors=imm_hmap_color,
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.6, 
          breaks=seq(0, 1, length.out=257))


##..Add cell type to each event..##
myeloid_clust<-c(1,2,3,11,12,13,14,21,22,23,24,25,30,31,32,33,34,35,36,41,42,44,45,51,52,53,54,55,56,57,61,62,63,64,65,71,72,73,74,75,81,82,83,84,91,92,93)
CD4_clust<-c(26,27,28,38,39,47,48,49,50,58,59,60,68,69,70,76,77,79,85,86,94,95,96)
CD8_clust<-c(6,15,16,29,40,66,78,80,87,88,89,90,97,98,99,100)
B_clust<-c(7,8,9,10,17,18,19,20)
treg_clust<-c(67)
neut_clust<-c(43,46)
other_clust<-c(4,37)
mast_clust=c(5)
imm_pheno<-fs_clusters_imm
imm_pheno<-replace(imm_pheno,imm_pheno %in% myeloid_clust,"myeloid")
imm_pheno<-replace(imm_pheno,imm_pheno %in% CD4_clust,"CD4_T")
imm_pheno<-replace(imm_pheno,imm_pheno %in% CD8_clust,"CD8_T")
imm_pheno<-replace(imm_pheno,imm_pheno %in% B_clust,"B_cell")
imm_pheno<-replace(imm_pheno,imm_pheno %in% treg_clust,"Treg")
imm_pheno<-replace(imm_pheno,imm_pheno %in% mast_clust,"mast")
imm_pheno<-replace(imm_pheno,imm_pheno %in% neut_clust,"neutrophil")
imm_pheno<-replace(imm_pheno,imm_pheno %in% other_clust,"imm_other")
data_immune$imm_pheno<-imm_pheno
data_immune[mast_rows,]$imm_pheno<-"mast"
data_immune[neut_rows,]$imm_pheno<-"neutrophil"

#manually annotate gdT cells
CD3_thresh<-mean(data_immune[data_immune$imm_pheno=="CD4_T",]$CD3)
gdT_row_idx<-which(data_immune$gdTCR > 0.5 & data_immune$CD3 >= CD3_thresh)
gdT_rows<-rownames(data_immune[gdT_row_idx,])
data_immune[gdT_rows,]$imm_pheno<-"gdT_cell"


                  ###################################################
                  ####..FlowSOM Cluster round 3: Myeloid Cells..#####
                  ###################################################

myeloid<-c("myeloid") #can append if wanting to include other pops
data_myeloid<-data_immune[data_immune$imm_pheno %in% myeloid,]
ff_myeloid <- flowFrame(exprs = data.matrix(data_myeloid[,-c(47,48,49,50)]), desc = list(FIL = 1))

clusterChannels_3=c("CD14","CD16","CD11c","CD11b","CD68","CD163","CD206","CD209")

##..Run FlowSOM random seed for reproducibility..##

set.seed(321)
out_fSOM_myeloid <- FlowSOM::ReadInput(ff_myeloid, transform = FALSE, scale = FALSE, compensate = FALSE)
out_fSOM_myeloid <- FlowSOM::BuildSOM(out_fSOM_myeloid, colsToUse = clusterChannels_3)
out_fSOM_myeloid <- FlowSOM::BuildMST(out_fSOM_myeloid)
labels_myeloid <- out_fSOM_myeloid$map$mapping[,1]

out_fSOM_myeloid <- UpdateNodeSize(out_fSOM_myeloid, reset = TRUE)

FlowSOM::PlotStars(out_fSOM_myeloid, view = "grid", markers = clusterChannels_3)
out_fSOM_myeloid <- UpdateNodeSize(out_fSOM_myeloid)
FlowSOM::PlotStars(out_fSOM_myeloid, view = "MST", markers = clusterChannels_3)
FlowSOM::PlotStars(out_fSOM_myeloid, view="tSNE",markers = clusterChannels_3)

##..Visualize initial FlowSOM clusters output on heatmap..##

# Get FlowSOM cluster assignments and append to matrix of percentile normalized data
fs_clusters_myeloid<-out_fSOM_myeloid[["map"]][["mapping"]]
fs_clusters_myeloid<-fs_clusters_myeloid[,1]
data_fs_clusters_myeloid <- cbind(data_myeloid, fs_clusters_myeloid)

# go through all clusters and calculate mean for every channel
hm_allclusters_myeloid <- matrix(, nrow = length(unique(fs_clusters_myeloid)), ncol = length(clusterChannels_3))
for(i in 1:length(unique(fs_clusters_myeloid))) {
  temp_mat <- data_fs_clusters_myeloid[data_fs_clusters_myeloid[,"fs_clusters_myeloid"] == i, clusterChannels_3]
  hm_allclusters_myeloid[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# add names to rows and cols
rownames(hm_allclusters_myeloid) <- paste("cluster", 1:length(unique(fs_clusters_myeloid)), sep = "")
colnames(hm_allclusters_myeloid) <- clusterChannels_3  
hm_allclusters_myeloid

# plot heatmap of all clusters

heatmap.2(hm_allclusters_myeloid, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          #col = colorRampPalette(rev(brewer.pal(11,"Spectral")))(100),
          col = viridis(256),
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.9, margins = c(8,14),
          breaks=seq(0, 1, length.out=257))


##..Meta-cluster..##

# try the suggested automatic metaclustering method for a hint for k
auto_meta_myeloid <- MetaClustering(out_fSOM_myeloid$map$codes, method = "metaClustering_consensus", max = 20)
max(auto_meta_myeloid)

# do a manual metaclustering
chosen_k=20
set.seed(321)
out_meta_myeloid <- FlowSOM::metaClustering_consensus(out_fSOM_myeloid$map$codes, k = chosen_k)
meta_results_myeloid <- out_meta_myeloid[labels_myeloid] 

##..Visualize FlowSOM metaclusters output on heatmap..##

# make combined expression matrix
data_meta_clusters_myeloid <- cbind(data_myeloid, meta_results_myeloid)

# check the amount of cells in each cluster
table(data_meta_clusters_myeloid$meta_results_myeloid)

# go through all clusters and calculate mean for every channel
hm_metaclusters_myeloid <- matrix(, nrow = chosen_k, ncol = length(clusterChannels_3))
for(i in 1:chosen_k) {
  temp_mat <- data_meta_clusters_myeloid[data_meta_clusters_myeloid[,"meta_results_myeloid"] == i, clusterChannels_3]
  hm_metaclusters_myeloid[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(hm_metaclusters_myeloid) <- paste("cluster", 1:chosen_k, sep = "")
colnames(hm_metaclusters_myeloid) <- clusterChannels_3 
hm_metaclusters_myeloid 

# make a metacluster heatmap
myeloid_hmap_color<-sort(unique(meta_results_myeloid))
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster1,"#F8766D")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster2,"#00C0B2")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster3,"#00BA4F")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster4,"#87AB01")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster5,"#D277FF")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster6,"#E7851E")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster7,"#29A3FF")
myeloid_hmap_color<-replace(myeloid_hmap_color,myeloid_hmap_color %in% cluster8,"#00BCD6")

heatmap.2(hm_metaclusters_myeloid, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both"),
          trace = "none",
          col = viridis(256),
          sepwidth=c(0.01,0.01),
          sepcolor="grey",
          colsep=0:ncol(hm_metaclusters_myeloid),
          rowsep=0:nrow(hm_metaclusters_myeloid),
          RowSideColors=myeloid_hmap_color,
          density.info = 'none',
          key.title = '',
          cexRow = 0.6, cexCol = 0.6, 
          breaks=seq(0, 1, length.out=257))

cluster1<-c(12) #CD14_Mono
cluster2<-c(1,8,9) #CD16_CD14_Mono
cluster3<-c(15,16,17) #CD68_Mac
cluster4<-c(11,13) #CD11c_DC/Mono
cluster5<-c(19,20) #CD209_DC
cluster6<-c(5,6,7,10,14) #CD11b/c_CD206_Mac/Mono
cluster7<-c(2,3,4) #CD163_Mac
cluster8<-c(18) #CD206_Mac


                  ###################################################
                  ##########..Phenograph: Myeloid Cells..###########
                  ###################################################

# run phenograph on cluster channels
library(Rphenograph)
phenograph<-Rphenograph(as.matrix(data_myeloid[,clusterChannels_3]))
pheno_clusters<-as.data.frame(factor(membership(phenograph[[2]])))
colnames(pheno_clusters)<-'phenograph'

# Visualize Phenograph clusters on heatmap and meta-cluster

# make combined expression matrix
data_pheno_myeloid<-cbind(data_myeloid,pheno_clusters)

# check the amount of cells in each cluster
table(data_pheno_myeloid$phenograph)

# go through all clusters and calculate mean for every channel
hm_phenoclusters_myeloid <- matrix(, nrow = 36, ncol = length(clusterChannels_3))
for(i in 1:36) {
  temp_mat <- data_pheno_myeloid[data_pheno_myeloid[,"phenograph"] == i, clusterChannels_3]
  hm_phenoclusters_myeloid[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}

# rename
rownames(hm_phenoclusters_myeloid) <- paste("cluster", 1:36, sep = "")
colnames(hm_phenoclusters_myeloid) <- clusterChannels_3 
hm_phenoclusters_myeloid 

heatmap.2(hm_phenoclusters_myeloid, 
          scale = "none",
          Colv = T, Rowv = T,
          hclustfun = hclust,
          dendrogram = c("both","row","column","none"),
          trace = "none",
          col = magma(256),
          density.info = 'none',
          key.title = '',
          lhei = c(1,7),
          cexRow = 0.6, cexCol = 0.9)
          #breaks=seq(0, 1, length.out=101))


myeloid_pheno<-meta_results_myeloid
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster1,"CD14_Mono")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster2,"CD16_CD14_Mono")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster3,"CD68_Mac")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster4,"CD11c_DC/Mono")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster5,"CD209_DC")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster6,"CD11b/c_CD206_Mac/Mono")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster7,"CD163_Mac")
myeloid_pheno<-replace(myeloid_pheno,myeloid_pheno %in% cluster8,"CD206_Mac")
data_myeloid$myeloid_pheno<-myeloid_pheno


              ##################################################################################
              #####..Concatenate into one matrix with terminal cell-type ID for all cells..#####
              ##################################################################################

data_gran_norm$cell_type<-data_gran_norm$lineage
data_gran_norm[rownames(data_immune),]$cell_type<-data_immune$imm_pheno
data_gran_norm[rownames(data_myeloid),]$cell_type<-data_myeloid$myeloid_pheno
data_gran$cell_type<-data_gran_norm$cell_type

#add details for major cell type (lymphocyte, myeloid, granulocyte, fibroblast, endothelial, epithelial)
myeloid<-c(unique(data_myeloid$myeloid_pheno))
lymphocyte<-c("CD8_T","CD4_T","B_cell","Treg","gdT_cell")
granulocyte<-c("mast","neutrophil")
imm_other<-c("imm_other")
nonimmune<-c("fibroblast","endothelial","epithelial")


data_gran_norm<-data_gran_norm %>% mutate(cell_lin=case_when(data_gran_norm$cell_type %in% myeloid ~ "myeloid",
                                                             data_gran_norm$cell_type %in% lymphocyte ~ "lymphocyte",
                                                             data_gran_norm$cell_type %in% granulocyte ~ "granulocyte",
                                                             data_gran_norm$cell_type %in% imm_other ~ "other",
                                                             data_gran_norm$cell_type %in% nonimmune ~ "nonimmune")) 
data_gran<-data_gran %>% mutate(cell_lin=case_when(data_gran$cell_type %in% myeloid ~ "myeloid",
                                                   data_gran$cell_type %in% lymphocyte ~ "lymphocyte",
                                                   data_gran$cell_type %in% granulocyte ~ "granulocyte",
                                                   data_gran$cell_type %in% imm_other ~ "other",
                                                   data_gran$cell_type %in% nonimmune ~ "nonimmune")) 

#add dummy/numerical indicators for easier indexing in matlab
data_gran<-data_gran %>% mutate(lintype_num=case_when(data_gran$cell_lin == "myeloid" ~ 1,
                                                   data_gran$cell_lin == "lymphocyte" ~ 2,
                                                   data_gran$cell_lin == "granulocyte" ~ 3,
                                                   data_gran$cell_lin == "other" ~ 4,
                                                   data_gran$cell_lin == "nonimmune" ~ 5))
data_gran_norm$lintype_num<-data_gran$lintype_num

write.csv(data_gran_norm, file="granA_cellpheno_CS-asinh-norm.csv",row.names = FALSE) #save the annotated asinh-99th percentile scaled data
write.csv(data_gran, file="granA_cellpheno_CS.csv",row.names = FALSE) #save the annotatetd CS, but untransformed data
