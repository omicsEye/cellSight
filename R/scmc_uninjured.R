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


pbmc.nonwound1.data <- Read10X(paste0("data/Nonwound1/filtered_feature_bc_matrix", sep= "")) 

# Create a Seurat object
pbmc.nonwound1 <- CreateSeuratObject(counts = pbmc.nonwound1.data, project = "3", min.cells = 3, min.features = 200)
pbmc.nonwound1$stim <- "Uninjured_1"
pbmc.nonwound1 <- subset(pbmc.nonwound1, subset = nFeature_RNA > 500)
pbmc.nonwound1 <- NormalizeData(pbmc.nonwound1, verbose = FALSE)
pbmc.nonwound1 <- FindVariableFeatures(pbmc.nonwound1, selection.method = "vst", nfeatures = 2000)  

pbmc.nonwound2.data <- Read10X(paste0("data/Nonwound2/filtered_feature_bc_matrix", sep= "")) 

# Create a Seurat object
pbmc.nonwound2 <- CreateSeuratObject(counts = pbmc.nonwound2.data, project = "4", min.cells = 3, min.features = 200)
pbmc.nonwound2$stim <- "Uninjured_2"
pbmc.nonwound2 <- subset(pbmc.nonwound2, subset = nFeature_RNA > 500)
pbmc.nonwound2 <- NormalizeData(pbmc.nonwound2, verbose = FALSE)
pbmc.nonwound2 <- FindVariableFeatures(pbmc.nonwound2, selection.method = "vst", nfeatures = 2000)  


object.list <- FindIntegrationAnchors(object.list = list(pbmc.nonwound1, pbmc.nonwound2), dims = 1:20)
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
ggsave("all_samples_uninjured.pdf", plot=dimension, width = 7.2, height = 4, units = "in", dpi = 350)

DefaultAssay(object.list) <- "RNA"
nk.markers <- FindConservedMarkers(object.list, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

Feature_dist<-FeaturePlot(object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                   "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                   "Fabp4","Ptprc","Pecam1","CD34"), min.cutoff = "q9",combine= FALSE)
plots <- lapply(X = Feature_dist, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
feature_distri <- CombinePlots(plots = Feature_dist,ncol =3) + theme(aspect.ratio = 1)
ggsave(paste0("feature_dist_uninjured.pdf", sep=""), plot=feature_distri,width = 15, height = 11)

Feature_dist_violin <- VlnPlot(object = object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                                 "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                                 "Fabp4","Ptprc","Pecam1","CD34"),combine =FALSE)
plots <- lapply(X = Feature_dist_violin, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
feature_distri <- CombinePlots(plots = Feature_dist_violin,ncol =3) + theme(aspect.ratio = 1)
ggsave(paste0("feature_dist_violoin_uninjured.pdf", sep=""), plot=feature_distri,width = 15, height = 11)
