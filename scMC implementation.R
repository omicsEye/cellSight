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
##This part is for the Rano
setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
#This part is for Dr.Rahnavard
#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

#data_dir <- '/Users/Rano/Desktop/Single_Cell_Wound/'
data_dir <- '~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data'

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
  
  pbmc.nonwound1.data <- Read10X(paste0("data/Nonwound1/filtered_feature_bc_matrix", sep= "")) 
  
  # Create a Seurat object
  pbmc.nonwound1 <- CreateSeuratObject(counts = pbmc.nonwound1.data, project = "3", min.cells = 3, min.features = 200)
  pbmc.nonwound1$stim <- "Uninjured_1"
  pbmc.nonwound1 <- subset(pbmc.nonwound1, subset = nFeature_RNA > 500)
  pbmc.nonwound1 <- NormalizeData(pbmc.nonwound1, verbose = FALSE)
  pbmc.nonwound1 <- FindVariableFeatures(pbmc.nonwound1, selection.method = "vst", nfeatures = 2000)  
  ### For testing
  #pbmc.nonwound1 <- pbmc
  
  
  pbmc.nonwound2.data <- Read10X(paste0("data/Nonwound2/filtered_feature_bc_matrix", sep= "")) 
  
  # Create a Seurat object
  pbmc.nonwound2 <- CreateSeuratObject(counts = pbmc.nonwound2.data, project = "4", min.cells = 3, min.features = 200)
  pbmc.nonwound2$stim <- "Uninjured_2"
  pbmc.nonwound2 <- subset(pbmc.nonwound2, subset = nFeature_RNA > 500)
  pbmc.nonwound2 <- NormalizeData(pbmc.nonwound2, verbose = FALSE)
  pbmc.nonwound2 <- FindVariableFeatures(pbmc.nonwound2, selection.method = "vst", nfeatures = 2000)  
  
  
  object.list <- FindIntegrationAnchors(object.list = list(pbmc.nonwound1,pbmc.nonwound2,pbmc.wound1,pbmc.wound2), dims = 1:20)
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
  dim_plot <- DimPlot(object.list, reduction = "umap", split.by = "stim")
  dim_plot$data$stim <- factor(x = dim_plot$data$stim, levels = c("Uninjured_1","Uninjured_2", "Injured_1","Injured_2")) # change the order of the factor levels
  dim_plot # print again the plot
  #setwd("C:/Users/ranoj/Desktop/Single_Cell_output/")
  ggsave("dim_plot.pdf", plot=dimension, width = 7.2, height = 4, units = "in", dpi = 350)
  dimension <- DimPlot(object.list, reduction = "umap", split.by = "stim")
  ggsave("all_samples.pdf", plot=dimension, width = 7.2, height = 4, units = "in", dpi = 350)
  
  
  ###transfer cell_types label using 1 dataset as reference and 1 as query
  skin.query <- pbmc.wound1
  skin.integrated <- pbmc
  
  
  skin.anchors <- FindTransferAnchors(reference = object.list, query = skin.query,
                                          dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = skin.anchors, refdata = object.list$Celltype,
                              dims = 1:30)
  skin.query <- AddMetaData(skin.query, metadata = predictions)
  skin.integrated <- RunUMAP(object.list, dims = 1:30, reduction = "pca", return.model = TRUE)
  skin.query <- MapQuery(anchorset = skin.anchors, reference = skin.integrated, query = skin.query,
                             refdata =list(celltype = "Celltype"), reference.reduction = "pca", reduction.model = "umap")
  skin.query.tsne <- MapQuery(anchorset = skin.anchors, reference = skin.integrated, query = skin.query,
                         refdata =list(celltype = "Celltype"), reference.reduction = "tsne", reduction.model = "tsne")
  p1 <- DimPlot(skin.integrated, reduction = "umap",group.by = "Celltype", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(skin.query, reduction = "ref.umap", group.by = "predicted.celltype",label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p1 + p2
  p3 <- DimPlot(skin.integrated, reduction = "ref.tsne",group.by = "Celltype", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  
  
  DefaultAssay(object.list) <- "RNA"
  nk.markers <- FindConservedMarkers(object.list, ident.1 = 0 , grouping.var = "stim", verbose = FALSE)
  head(nk.markers)
 
  Feature_dist<-FeaturePlot(object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                            "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                            "Fabp4","Ptprc","Pecam1","CD34"), min.cutoff = "q9",combine= FALSE)
   plots <- lapply(X = Feature_dist, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
   feature_distri <- CombinePlots(plots = Feature_dist,ncol =3) + theme(aspect.ratio = 1)
   ggsave(paste0("feature_dist.pdf", sep=""), plot=feature_distri,width = 15, height = 11)
   
   Feature_dist_violin <- VlnPlot(object = object.list, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                      "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                      "Fabp4","Ptprc","Pecam1","CD34"),combine =FALSE)
   plots <- lapply(X = Feature_dist_violin, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
   feature_distri <- CombinePlots(plots = Feature_dist_violin,ncol =3) + theme(aspect.ratio = 1)
   ggsave(paste0("feature_dist_violoin.pdf", sep=""), plot=feature_distri,width = 15, height = 11)
   #Finding the marker genes
  
object.list <- list("pbmc.wound1","pbmc.wound2")
  data.input <- object.list # a list of count data matrix, one per dataset
  sample.name <- names(data.input)
  object.list <- vector("list", length(sample.name))
  names(object.list) <- sample.name
  pbmc.combined <- merge(pbmc.wound1, y = c(pbmc.wound2, pbmc.nonwound1,pbmc.nonwound2), add.cell.ids = c("1", "2", "3","4"), project = "Mouse")
  #pbmc.combined
  #setwd("/Users/ranoj/Desktop/Single_Cell_output/")
  save(pbmc.combined, file = "pbmc_combined.RData")
  
  
  load("pbmc.combined.RData")  
  
  data.input <- pbmc.combined # a list of count data matrix, one per dataset
  sample.name <- names(data.input)
  
  object.list <- SplitObject(pbmc.combined, split.by = "orig.ident")
  object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  
  features <- SelectIntegrationFeatures(object.list = object.list)
  object.list <- lapply(X = object.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  features.integration = identifyIntegrationFeatures(object.list)
  object.list <- identifyConfidentCells(object.list, features.integration)
  combined <- RunscMC(pbmc.combined)  
  
  ###Running for 1 injured and 1 uninjured.
  
  
  # object.list <- FindIntegrationAnchors(object.list = list(pbmc,pbmc_nonwound1), dims = 1:20)
  # object.list <- IntegrateData(anchorset = object.list, dims = 1:20)
  # DefaultAssay(object.list) <- "integrated"
  # 
  # object.list <- ScaleData(object.list, verbose = FALSE)
  # object.list <- RunPCA(object.list, npcs = 30, verbose = FALSE)
  # object.list <- RunUMAP(object.list, reduction = "pca", dims = 1:30, verbose = FALSE)
  # object.list <- RunTSNE(object.list, reduction = "pca", dims = 1:30, verbose = FALSE)
  # p1 <- DimPlot(object.list, reduction = "umap", group.by = "tech")
  # p2 <- DimPlot(object.list, reduction = "umap", group.by = "Celltype", label = TRUE, repel = TRUE) +
  #   NoLegend()
  # p1 + p2  
  # p3 <- DimPlot(object.list, reduction = "umap", group.by = "tech")
  # p4 <- DimPlot(object.list, reduction = "tsne", group.by = "Celltype", label = TRUE, repel = TRUE) +
  #   NoLegend()
  