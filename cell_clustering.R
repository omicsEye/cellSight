install.packages('Seurat')
install.packages('dplyr')
install.packages('Matrix')
install.packages('ggpubr')
install.packages('ggplots')
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)

setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

data_dir <- 'C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data'
#data_dir <- '~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data'

list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


# loop over samples and save figures and results 
#sample_list <- c("Wound1", "Wound2")

pbmc.data <- Read10X(data.dir = paste0("data/Nonwound1/filtered_feature_bc_matrix", sep= "")) 
pbmc.data <-Read10X(data.dir = "data/Nonwound1/filtered_feature_bc_matrix")
pbmc.data <- Read10X(data.dir = paste0("data/Wound2/filtered_feature_bc_matrix", sep= ""))  
pbmc.data <- Read10X(data.dir = paste0("data/Nonwound2/filtered_feature_bc_matrix", sep= ""))
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,
                               thresh.use = 0.25)
library(future)
library(scCATCH)
clu_markers <- findmarkergene(object  = pbmc,
                                species = 'Mouse',
                                cluster = 'All',
                                marker = 'cellmatch',
                                cancer = NULL,
                                tissue = NULL,
                                cell_min_pct = 0.25,
                                logfc = 0.25,
                                pvalue = 0.05)
#obj <- findmarkergene(object = pbmc, species = "Mouse", marker = cellmatch, tissue = "Blood")
clu_markers <- findmarkergenes(dif.PAC, species = "Human", cluster = 'All', match_CellMatch = TRUE, cancer = "Pancreatic Ductal Adenocarcinomas", tissue = "Pancreas", cell_min_pct = 0.25, logfc = 0.25, pvalue = 0.05)
clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Mouse',
                   cancer = NULL,
                   tissue = 'Adipose tissue')
new.cluster.ids <- c("N/A","Macrophage", "B-cell", "Neurons", "N/A","Adipocytes" )
new.cluster.ids <- c("Beta","Fibroblast","B-cell","Neurons","Trigeminal neurons","Endothelial cells",
                     "Smooth Muscle cells","Keratinocytes","Keratinocytes","Macrophages","Smooth Muscle cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc@meta.data <- pbmc@meta.data %>%
  mutate(Cell_type = case_when(
    seurat_clusters == 0 ~ "N/A",
    seurat_clusters == 1 ~ "Macrophages",
    seurat_clusters == 2 ~ "B-cell",
    seurat_clusters == 3 ~ "Neurons",
    seurat_clusters == 4 ~ "N/A",
    seurat_clusters == 5 ~ "Adipocytes"
    
  ))
test<-t(pbmc.data)

write.table(pbmc@meta.data,"meta_data_wound2.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(test,"data_wound2.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)

pbmc@meta.data <- pbmc@meta.data %>%
  mutate(Cell_type = case_when(
    seurat_clusters == 0 ~ "Beta",
    seurat_clusters == 1 ~ "Fibroblast",
    seurat_clusters == 2 ~ "B-cell",
    seurat_clusters == 3 ~ "Neurons",
    seurat_clusters == 4 ~ "Trigeminal neurons",
    seurat_clusters == 5 ~ "Endothelial cells",
    seurat_clusters == 6 ~ "Smooth Muscle cells",
    seurat_clusters == 7 ~ "Keratinocytes",
    seurat_clusters == 8 ~ "Keratinocytes",
    seurat_clusters == 9 ~ "Macrophages",
    seurat_clusters == 10 ~ "Smooth Muscle cells",
    
  ))

test<-t(pbmc.data)

write.table(pbmc@meta.data,"meta_data_wound2.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(test,"data_wound2.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)


wound1 <- pbmc
wound1.data <- pbmc.data 

wound2 <- pbmc
wound2.data <- pbmc.data 

unwound1 <- pbmc
unwound1.data <- pbmc.data
# Sample_1$Sample <- "Sample_1"
}
