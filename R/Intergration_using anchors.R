# script to integrate scRNA-Seq datasets to correct for batch effects
setwd("C:/Users/ranoj/Desktop/Single_Cell_output/")


# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
BiocManager::install("SingleR",force = T)
library(SingleR)
##Removing all objects from the memory
rm(list = ls())
# get data location
dirs <- list.dirs(path = 'data/', recursive = F, full.names = F)

x<- 'Nonwound1_skin_filtered_feature_bc_matrix'
for(x in dirs){
  name <- gsub('_filtered_feature_bc_matrix','', x)

  cts <- Read10X(paste0('data/',x, sep= ""))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}




# merge datasets

merged_seurat <- merge(Nonwound1_uninjuredskin, y = c(Nonwound2_uninjuredskin, Wound1_injuredskin, Wound2_injuredskin),
                       add.cell.ids = ls()[4:7],
                       project = 'Injured_mouse')


merged_seurat

# QC & filtering -----------------------

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

#


# filtering

merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 5)

merged_seurat_filtered

merged_seurat




# perform standard workflow steps to figure out if we see any batch effects --------
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)

df1 <- seurat.integrated@meta.data[seurat.integrated@meta.data$Patient == 'Nonwound1', ]
df2 <- seurat.integrated@meta.data[seurat.integrated@meta.data$Patient == 'Nonwound2', ]
df3 <- seurat.integrated@meta.data[seurat.integrated@meta.data$Patient == 'Wound1', ]
df4 <- seurat.integrated@meta.data[seurat.integrated@meta.data$Patient == 'Wound2', ]


write.table(df1, file='Nonwound_1_matched.tsv', quote=FALSE, sep='\t')
write.table(df2, file='Nonwound_2_matched.tsv', quote=FALSE, sep='\t')
write.table(df3, file='Wound_1_matched.tsv', quote=FALSE, sep='\t')
write.table(df4, file='Wound_2_matched.tsv', quote=FALSE, sep='\t')

# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))


grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
saveRDS(seurat.integrated, file = "seurat_integrated.rds")

set.seed(1234)

install.packages("harmony")
library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)





# load dataset
skin.filtered <- merged_seurat


# QC and filtering


# standard workflow steps
skin.filtered <- NormalizeData(skin.filtered)
skin.filtered <- FindVariableFeatures(skin.filtered)
skin.filtered <- ScaleData(skin.filtered)
skin.filtered <- RunPCA(skin.filtered)
ElbowPlot(skin.filtered)
skin.filtered.umap <- RunUMAP(skin.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(skin.filtered.umap, reduction = 'umap', group.by = 'Type')
before.simple <- DimPlot(skin.filtered.umap, reduction = 'umap')

p1 <- DimPlot(skin.filtered.umap, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(skin.filtered.umap, reduction = 'umap', group.by = 'Type',
              cols = c('red','green','blue'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)

# run Harmony -----------
skin.harmony <- skin.filtered %>%
  RunHarmony(group.by.vars = 'Type', plot_convergence = FALSE)

skin.harmony@reductions

skin.harmony.embed <- Embeddings(skin.harmony, "harmony")
skin.harmony.embed[1:50,1:50]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
skin.harmony <- skin.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.5)

# visualize 
after <- DimPlot(skin.harmony, reduction = 'umap', group.by = 'Type')
after.simple <- DimPlot(skin.harmony, reduction = 'umap')


before|after
before.simple|after.simple

clusters <- DimPlot(skin.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(skin.harmony, reduction = 'umap', group.by = 'Type')

###Reordering for plotting
skin.harmony$Type <- factor(x = skin.harmony$Type, levels = c("uninjuredskin", "injuredskin"))


groups <- DimPlot(skin.harmony, reduction = 'umap', group.by = 'seurat_clusters',split.by = 'Type', label = TRUE)

condition|clusters

skin.harmony.markers <-FindAllMarkers(skin.harmony,only.pos = FALSE, min.pct = 0.01)
skin.harmony.markers<- skin.harmony.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
folder <- "markers"

if (file.exists(folder)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(folder)
  
}
write.table(skin.harmony.markers, file=paste0("markers/harmony_all.tsv"), quote=FALSE, sep='\t', col.names = NA)
View(skin.harmony@meta.data)
Feature_dist <- VlnPlot(object = skin.harmony, features =c("Ly6G","Adgre1","Ptprc","Pdgfra","Pdgfrb","Col1a1",
                                                   "Krt14","Krt10","Krt5","Plin1","Adipoq","Pparg",
                                                   "Fabp4","Ptprc","Pecam1","CD25"),combine =FALSE)
plots <- lapply(X = Feature_dist, FUN = function(x) x + theme(axis.text = element_text(size = 4)))
feature_distri <- CombinePlots(plots = Feature_dist,ncol =3) + theme(aspect.ratio = 1)
ggsave(paste0("feature_dist_.pdf", sep=""), plot=feature_distri,width = 15, height = 11)
#Finding the marker genes

ref.data <-celldex::ImmGenData()

BiocManager::install("SingleCellExperiment")
BiocManager::install('scran')

##convert the seyurat object to singlecell object 
sce <- as.SingleCellExperiment(DietSeurat(skin.harmony))
sce

ref.data.main <- SingleR(test = sce,assay.type.test = 1,ref = ref.data,labels = ref.data$label.main,de.method = "wilcox")
#The second reference is called fine, which has finer references to the celltypes
ref.data.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref.data,labels = ref.data$label.fine)
table(ref.data.main$pruned.labels)
table(ref.data.fine$pruned.labels)
skin.harmony@meta.data$ref.data.main <- ref.data.main$pruned.labels
skin.harmony@meta.data$ref.data.fine <- ref.data.fine$pruned.labels

df_fibro<-  df_full[df_full$ref.data.main== 'Fibroblasts', ]
