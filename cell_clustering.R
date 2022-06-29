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
#install.packages(pkgs = 'devtools')
#devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)
##This part is for the Rano
setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
#This part is for Dr.Rahnavard
#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

data_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/"
#data_dir <- '~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data'

list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


# loop over samples and save figures and results 
#sample_list <- c("Wound1", "Wound2", "Nonwound1", "Nonwound2")
sample <- "Nonwound2"
for (sample in sample_list){
  
  print(sample)
  # read data  
  pbmc.data <- Read10X(paste0("data/", sample, "/filtered_feature_bc_matrix", sep= "")) 
  
  # Create a Seurat object
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  
  # Shows the mitrochrondial percentage in the data
  
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  #output to a different path
  setwd("C:/Users/ranoj/Desktop/Single_Cell_output/")
  QC_VlnPlot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(paste0("analysis/figures/QC_Plots/QC_VlnPlot_", sample,".pdf", sep=""), plot=QC_VlnPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(paste0("analysis/figures/QC_Plots/QC_Scatter_", sample,".pdf", sep=""), plot=plot1 + plot2, width = 7.2, height = 4, units = "in", dpi = 350)
  
  # Normalization of the data
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- NormalizeData(pbmc)
  
  # Feature Selection
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(pbmc), 10)
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(pbmc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(paste0("analysis/figures/QC_Plots/Top_10Genes_", sample,".pdf", sep=""), plot=plot1 + plot2, width = 14, height = 4, units = "in", dpi = 350)
  
  # 
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  
  # make matrix and do a heatmap use phetamap for loading and scores in PCA for 5 PCs and 20 Features
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
  PCA_score_1 <- VizDimLoadings(pbmc, dims = 1, reduction = "pca") +  theme(axis.text = element_text(size = 6))
  PCA_score_2 <- VizDimLoadings(pbmc, dims = 2, reduction = "pca")+  theme(axis.text = element_text(size = 6))
  ggsave(paste0("analysis/figures/QC_Plots/PCA_Scores_", sample,".pdf", sep=""), plot=PCA_score_1|PCA_score_2, width = 7.2, height = 4, units = "in", dpi = 350)
  
  #PCA plots
  QC_DimPlot <- DimPlot(pbmc, reduction = "pca")
  QC_DimHeatmap <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
  QC_DimHeatmap15 <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
  #Plot the Jackstraw score
  pbmc <- JackStraw(pbmc, num.replicate = 100)
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
  QC_JackStrawPlot <- JackStrawPlot(pbmc, dims = 1:15)
  ggsave(paste0("analysis/figures/QC_Plots/Elbow_plot", sample,".pdf", sep=""), plot=QC_JackStrawPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  #Elbow plot
  QC_ElbowPlot <- ElbowPlot(pbmc)
  ggsave(paste0("analysis/figures/QC_Plots/Elbow_plot", sample,".pdf", sep=""), plot=QC_ElbowPlot, width = 7.2, height = 4, units = "in", dpi = 350)
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = F)
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5)
  DimPlot(pbmc, label = T)
  head(Idents(pbmc), 5)
  #UMAP plot with clusters

  
  
  
  ## Load the Immogen dataset using celldex

  ref.data <-celldex::ImmGenData()
  
  BiocManager::install("SingleCellExperiment")
  BiocManager::install('scran')
  
  ##convert the seyurat object to singlecell object 
  sce <- as.SingleCellExperiment(DietSeurat(pbmc))
  sce
  #Immogen data set has 2 ref database one being main, which has all the main celltypes 
  ref.data.main <- SingleR(test = sce,assay.type.test = 1,ref = ref.data,labels = ref.data$label.main,de.method = "wilcox")
  #The second reference is called fine, which has finer references to the celltypes
  ref.data.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref.data,labels = ref.data$label.fine)
  table(ref.data.main$pruned.labels)
  table(ref.data.fine$pruned.labels)
  pbmc@meta.data$ref.data.main <- ref.data.main$pruned.labels
  pbmc@meta.data$ref.data.fine <- ref.data.fine$pruned.labels
  
  meta_nonwound1 <- pbmc@meta.data
  colmn <- paste("col", 1:3)
  meta_nonwound1<- meta_nonwound1 %>% separate(ref.data.fine, sep = " ", into = colmn, remove = FALSE)
  meta_nonwound1 <- meta_nonwound1[ , -which(names(meta_nonwound1) %in% c("col 2","col 3"))]
  names(meta_nonwound1)[names(meta_nonwound1) == 'col 1'] <- "Cell_types"
  result <- meta_nonwound1 %>% 
    add_count(seurat_clusters, Cell_types,sort = T) %>%
    group_by(seurat_clusters) %>%
    mutate(Majority = Cell_types[n == max(n)][1]) %>%
    # do not keep temp var
    select(-n)
  labels <- pbmc@meta.data
  
  new.cluster.ids <- c("Fibroblasts", "Fibroblasts", "Fibroblasts", "Macrophages", "Fibroblasts", "Fibroblasts",
                       "Endothelial", "Stromal", "Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts","Fibroblasts")
  names(new.cluster.ids) <- levels(pbmc)
  pbmc <- RenameIdents(pbmc, new.cluster.ids)
  DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
 ###Remove this when done
  pbmc_nonwound2 <- pbmc
  pbmc_nonwound1 <- pbmc
  pbmc$Celltype <- Idents(pbmc)
  pbmc_nonwound1$Celltype <- Idents(pbmc_nonwound1)
  write.table(labels, file=paste0("analysis/data/Label_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  # pbmc <- SetIdent(pbmc, value = "ref.data.main")
  # DimPlot(pbmc, reduction = "umap" ,label = T , repel = T, label.size = 3) + NoLegend()
  # pbmc <- SetIdent(pbmc, value = "ref.data.fine")
  # pbmc <- SetIdent(pbmc, value = "ref.data.main")
  #cluster <-DimPlot(pbmc, reduction = "umap" ,label = T , repel = T, label.size = 3) + NoLegend()
  
  
  # hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
  # hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
  # dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
  # dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
  ##Shows the marker for cluster 2 only
  cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
  #pbmc@meta.data$markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
  head(cluster2.markers, n = 5)
  # displays the marker for every clusters compaared to all remaining cell
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  write.table(pbmc.markers, file=paste0("analysis/data/marker_gene_unique_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  Feature_dist <- VlnPlot(object = pbmc, features =c("CD68","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                     "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                     "Fabp4","Ptprc","Pecam1","CD34"),combine =FALSE)
  plots <- lapply(X = Feature_dist, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
  feature_distri <- CombinePlots(plots = Feature_dist,ncol =3) + theme(aspect.ratio = 1)
  ggsave(paste0("analysis/figures/QC_Plots/feature_dist_", sample,".pdf", sep=""), plot=feature_distri,width = 15, height = 11)
  #Finding the marker genes
  pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,
                                 thresh.use = 0.25)
  cluster_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +ylim(-10,15) +xlim(-10,15)
  ggsave(paste0("analysis/figures/QC_Plots/Cluster_without_label", sample,".pdf", sep=""), plot=cluster_plot, width = 7.2, height = 4, units = "in", dpi = 350)
  print(sample)
  #Saving the marker gene files for better understanding the process
  test <- pbmc.markers
  write.table(test, file=paste0("analysis/data/marker_gene_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  test<-as.matrix(pbmc.data)
  test <- t(test)
  final <- rownames(labels)
  final_data <- subset(test, rownames(test)%in% final)
  write.table(final_data, file=paste0("analysis/data/data_labelled_", sample,".tsv"), quote=FALSE, sep='\t', col.names = NA)
  # library(future)
  # library(scCATCH)
  # # clu_markers <- findmarkergene(object  = pbmc,
  # #                                 species = 'Mouse',
  # #                                 cluster = 'All',
  # #                                 marker = 'cellmatch',
  # #                                 cancer = NULL,
  # #                                 tissue = NULL,
  # #                                 cell_min_pct = 0.25,
  # #                                 logfc = 0.25,
  # #                                 pvalue = 0.05)
  # #obj <- findmarkergene(object = pbmc, species = "Mouse", marker = cellmatch, tissue = "Blood")
  # # clu_markers <- findmarkergenes(dif.PAC, species = "Human", cluster = 'All', match_CellMatch = TRUE, cancer = "Pancreatic Ductal Adenocarcinomas", tissue = "Pancreas", cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05)
  # # clu_ann <- scCATCH(object = clu_markers$clu_markers,
  # #                    species = 'Mouse',
  # #                    cancer = NULL,
  # #                    tissue = 'Adipose tissue')
  # if (sample == "Wound1" | sample == "Wound2"){
  # new.cluster.ids <- c("N/A","Macrophage", "B-cell", "Neurons", "N/A","Adipocytes" )
  # }
  # if (sample == "Nonwound1"){
  # new.cluster.ids <- c("Beta","Fibroblast","B-cell","Neurons","Trigeminal neurons","Endothelial cells",
  #                       "Smooth Muscle cells","Keratinocytes","Keratinocytes","Macrophages","Smooth Muscle cells")
  # }
  # if (sample == "Nonwound2"){
  #   new.cluster.ids <- c("Beta","Fibroblast","B-cell","Neurons","Trigeminal neurons","Endothelial cells",
  #                        "Smooth Muscle cells","Keratinocytes","Keratinocytes","Macrophages","Smooth Muscle cells","N/A","N/A")
  # }
  # names(new.cluster.ids) <- levels(pbmc)
  # pbmc <- RenameIdents(pbmc, new.cluster.ids)
  # cluster_plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+ylim(-8,15) +xlim(-10,10)
  # ggsave(paste0("analysis/figures/QC_Plots/Cluster_with_label", sample,".pdf", sep=""), plot=cluster_plot, width = 7.2, height = 4, units = "in", dpi = 350)
  # if (sample == "Wound1" | sample == "Wound2"){
  # pbmc@meta.data <- pbmc@meta.data %>%
  #   mutate(Cell_type = case_when(
  #     seurat_clusters == 0 ~ "N/A",
  #     seurat_clusters == 1 ~ "Macrophages",
  #     seurat_clusters == 2 ~ "B-cell",
  #     seurat_clusters == 3 ~ "Neurons",
  #     seurat_clusters == 4 ~ "N/A",
  #     seurat_clusters == 5 ~ "Adipocytes"
  # 
  #   ))
  # }
  # else if (sample == "Nonwound1") {
  #   pbmc@meta.data <- pbmc@meta.data %>%
  #     mutate(Cell_type = case_when(
  #       seurat_clusters == 0 ~ "Beta",
  #       seurat_clusters == 1 ~ "Fibroblast",
  #       seurat_clusters == 2 ~ "B-cell",
  #       seurat_clusters == 3 ~ "Neurons",
  #       seurat_clusters == 4 ~ "Trigeminal neurons",
  #       seurat_clusters == 5 ~ "Endothelial cells",
  #       seurat_clusters == 6 ~ "Smooth Muscle cells",
  #       seurat_clusters == 7 ~ "Keratinocytes",
  #       seurat_clusters == 8 ~ "Keratinocytes",
  #       seurat_clusters == 9 ~ "Macrophages",
  #       seurat_clusters == 10 ~ "Smooth Muscle cells"
  #     ))
  # }
#   else if (sample == "Nonwound2") {
#     pbmc@meta.data <- pbmc@meta.data %>%
#       mutate(Cell_type = case_when(
#         seurat_clusters == 0 ~ "Beta",
#         seurat_clusters == 1 ~ "Fibroblast",
#         seurat_clusters == 2 ~ "B-cell",
#         seurat_clusters == 3 ~ "Neurons",
#         seurat_clusters == 4 ~ "Trigeminal neurons",
#         seurat_clusters == 5 ~ "Endothelial cells",
#         seurat_clusters == 6 ~ "Smooth Muscle cells",
#         seurat_clusters == 7 ~ "Keratinocytes",
#         seurat_clusters == 8 ~ "Keratinocytes",
#         seurat_clusters == 9 ~ "Macrophages",
#         seurat_clusters == 10 ~ "Smooth Muscle cells",
#         seurat_clusters == 11 ~ "N/A",
#         seurat_clusters == 12 ~ "N/A"
#       ))
#   }
#   #Saving the data in the desired format
# 
#   test <- as.matrix(pbmc.data)
#   test<-t(test)
#   # Saving the data and metadata for tweedeverse analysis
#   write.table(pbmc@meta.data,paste0("analysis/data/meta_data_", sample, ".tsv", sep=""),  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
#   write.table(test,paste0("analysis/data/data_", sample, ".tsv", sep=""),  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
}


