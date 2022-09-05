install.packages('Seurat')
install.packages('dplyr')
install.packages('Matrix')
install.packages('ggpubr')
install.packages('ggplots')
#install.packages('SingleR')
BiocManager::install("celldex",force = T)
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install('CHETAH')

BiocManager::install("SingleR",force = T)
devtools::install_github("amsszlh/scMC")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(CHETAH)
library(celldex)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(cowplot)

install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')
library(multtest)
library(metap)

library(scMC)
#install.packages(pkgs = 'devtools')
#devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)
##This part is for Windows
setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")

##This is for MAC
##setwd("/Users/Rano/Desktop/Single_Cell_Wound/")
#This part is for Dr.Rahnavard
#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

data_dir <- '/Users/Rano/Desktop/Single_Cell_Wound/'
#data_dir <- '~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data'

list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


# loop over samples and save figures and results 
#sample_list <- c("Wound1", "Wound2", "Nonwound1", "Nonwound2")
#sample <- "Wound1"
#for (sample in sample_list){

#print(sample)
# read data  
pbmc.wound1.data <- Read10X(paste0("data/Wound1/filtered_feature_bc_matrix", sep= "")) 

# Create a Seurat object
pbmc.wound1 <- CreateSeuratObject(counts = pbmc.wound1.data, project = "1", min.cells = 3, min.features = 200)
pbmc.wound1$stim <- "Injured_1"
pbmc.wound1 <- subset(pbmc.wound1, subset = nFeature_RNA > 500)
pbmc.wound1 <- NormalizeData(pbmc.wound1, verbose = FALSE)
pbmc.wound1 <- FindVariableFeatures(pbmc.wound1, selection.method = "vst", nfeatures = 2000)


pbmc.wound2.data <- Read10X(paste0("data/Wound2/filtered_feature_bc_matrix", sep= "")) 

# Create a Seurat object
pbmc.wound2 <- CreateSeuratObject(counts = pbmc.wound2.data, project = "2", min.cells = 3, min.features = 200)
pbmc.wound2$stim <- "Injured_2"
pbmc.wound2 <- subset(pbmc.wound2, subset = nFeature_RNA > 500)
pbmc.wound2 <- NormalizeData(pbmc.wound2, verbose = FALSE)
pbmc.wound2 <- FindVariableFeatures(pbmc.wound2, selection.method = "vst", nfeatures = 2000)


object.list <- FindIntegrationAnchors(object.list = list(pbmc.wound1,pbmc.wound2), dims = 1:20)
object.list <- IntegrateData(anchorset = object.list, dims = 1:20)
DefaultAssay(object.list) <- "integrated"


# Run the standard workflow for visualization and clustering
object.list <- ScaleData(object.list, verbose = FALSE)
object.list <- RunPCA(object.list, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
object.list <- RunUMAP(object.list, reduction = "pca", dims = 1:20)
object.list <- FindNeighbors(object.list, reduction = "pca", dims = 1:20)
object.list <- FindClusters(object.list, resolution = 0.5)
p1 <- DimPlot(object.list, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object.list, reduction = "umap", label = TRUE)
plot_grid(p2, p1)
DimPlot(object.list, reduction = "umap", group.by = c("stim","ident"))

dimension <- DimPlot(object.list, reduction = "umap", split.by = "stim")
ggsave("wound_samples.pdf", plot=dimension, width = 7.2, height = 4, units = "in", dpi = 350)

DefaultAssay(object.list) <- "RNA"
nk.markers <- FindConservedMarkers(object.list, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

Feature_dist<-FeaturePlot(object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                   "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                   "Fabp4","Ptprc","Pecam1","CD34"), min.cutoff = "q9",combine= FALSE)
plots <- lapply(X = Feature_dist, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
feature_distri <- CombinePlots(plots = Feature_dist,ncol =3) + theme(aspect.ratio = 1)
ggsave(paste0("feature_dist_wound.pdf", sep=""), plot=feature_distri,width = 15, height = 11)

Feature_dist_violin <- VlnPlot(object = object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                                 "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                                 "Fabp4","Ptprc","Pecam1","CD34"),combine =FALSE)
plots <- lapply(X = Feature_dist_violin, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
feature_distri <- CombinePlots(plots = Feature_dist_violin,ncol =3) + theme(aspect.ratio = 1)
ggsave(paste0("feature_dist_violoin_wound.pdf", sep=""), plot=feature_distri,width = 15, height = 11)
#Finding the marker genes
