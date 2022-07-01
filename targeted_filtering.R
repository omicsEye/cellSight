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
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

BiocManager::install("SingleR",force = T)
library(monocle3)
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



setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")
sample <- "Nonwound2"

for (sample in sample_list)
  {
  print(sample)
  # read data  
  pbmc.data <- Read10X(paste0("data/", sample, "/filtered_feature_bc_matrix", sep= "")) 
  
  # Create a Seurat object
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  
  # Shows the mitrochrondial percentage in the data
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  #output to a different path
  setwd("~/Desktop/Single_Cell_output/")
  QC_VlnPlot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #ggsave(paste0("analysis/figures/QC_Plots/QC_VlnPlot_", sample,".pdf", sep=""), plot=QC_VlnPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #ggsave(paste0("analysis/figures/QC_Plots/QC_Scatter_", sample,".pdf", sep=""), plot=plot1 + plot2, width = 7.2, height = 4, units = "in", dpi = 350)
  
  # Normalization of the data
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- NormalizeData(pbmc)
  
  sample <- "Pdgfra"
  
  Pdgfra_expression = GetAssayData(object = pbmc, assay = "RNA", slot = "data")["Pdgfra",]
  pos_ids = names(which(Pdgfra_expression>0))
  
  pos_cells = subset(pbmc,cells=pos_ids)
  pos_cells <- FindVariableFeatures(pos_cells, selection.method = "vst", nfeatures = 2000)
  
  top10 <- head(VariableFeatures(pos_cells), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(pos_cells)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(paste0("analysis/figures/QC_Plots/Top_10Genes_", sample,".pdf", sep=""), plot=plot1 + plot2, width = 14, height = 4, units = "in", dpi = 350)
  
  # 
  all.genes <- rownames(pos_cells)
  pos_cells <- ScaleData(pos_cells, features = all.genes)
  
  pos_cells <- RunPCA(pos_cells, features = VariableFeatures(object = pos_cells))
  
  # make matrix and do a heatmap use phetamap for loading and scores in PCA for 5 PCs and 20 Features
  print(pos_cells[["pca"]], dims = 1:5, nfeatures = 5)
  PCA_score_1 <- VizDimLoadings(pos_cells, dims = 1, reduction = "pca") +  theme(axis.text = element_text(size = 6))
  PCA_score_2 <- VizDimLoadings(pos_cells, dims = 2, reduction = "pca")+  theme(axis.text = element_text(size = 6))
  ggsave(paste0("analysis/figures/QC_Plots/PCA_Scores_", sample,".pdf", sep=""), plot=PCA_score_1|PCA_score_2, width = 7.2, height = 4, units = "in", dpi = 350)
  
  #PCA plots
  QC_DimPlot <- DimPlot(pos_cells, reduction = "pca")
  QC_DimHeatmap <- DimHeatmap(pos_cells, dims = 1, cells = 500, balanced = TRUE)
  QC_DimHeatmap15 <- DimHeatmap(pos_cells, dims = 1:15, cells = 500, balanced = TRUE)
  #Plot the Jackstraw score
  pos_cells <- JackStraw(pos_cells, num.replicate = 100)
  pos_cells <- ScoreJackStraw(pos_cells, dims = 1:20)
  QC_JackStrawPlot <- JackStrawPlot(pos_cells, dims = 1:15)
  ggsave(paste0("analysis/figures/QC_Plots/Elbow_plot", sample,".pdf", sep=""), plot=QC_JackStrawPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  #Elbow plot
  QC_ElbowPlot <- ElbowPlot(pos_cells)
  ggsave(paste0("analysis/figures/QC_Plots/Elbow_plot", sample,".pdf", sep=""), plot=QC_ElbowPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  pos_cells <- FindNeighbors(pos_cells, dims = 1:10)
  pos_cells <- FindClusters(pos_cells, resolution = 0.5)
  pos_cells <- RunUMAP(pos_cells, dims = 1:10, verbose = F)
  pos_cells <- FindNeighbors(pos_cells, dims = 1:10)
  pos_cells <- FindClusters(pos_cells, resolution = 0.5)
  DimPlot(pos_cells, label = T)
  pos_cells.markers <- FindAllMarkers(pos_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top_2 <- pos_cells.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  write.table(pos_cells.markers, file=paste0("data_labelled_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  write.table(pos_cells.markers, file=paste0("data_top_2_labelled_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  head(Idents(pos_cells), 5)
  features <- c("CD68", "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1","Krt14", "Krt10", "Krt5","Plin1", "Adipoq", "Pparg", "Fabp4","Ptprc","Pecam1", "CD34")
  ridge <- RidgePlot(pos_cells, features = features, ncol = 3)
  ggsave(paste0("analysis/figures/ridge", sample,".pdf", sep="") ,plot=ridge, width = 7.2, height = 4, units = "in", dpi = 350)
  violin <- VlnPlot(pos_cells, features = features)
  feature_plot <- FeaturePlot(pos_cells, features = features)
  }