# load libraries ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

box_dir <- "~/Box/snRNA_CellRanger_Wound_nonWound/data/"
dirs <- list.dirs(path = box_dir,
                  recursive = F,
                  full.names = F)

obj_list <- dirs %>%
  set_names() %>%
  map(
    .f = function(x) {
      paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
        Read10X() %>%
        CreateSeuratObject()
    }
  )

obj_list$Nonwound1 %>% Cells() %>% length()
obj_list$Nonwound2 %>% Cells() %>% length()
obj_list$Wound1 %>% Cells() %>% length()
obj_list$Wound1 %>% Cells() %>% length()

# nonwound 1 qc -----------------------------------------------------------

obj_list$Nonwound1[["percent.mt"]] <- obj_list$Nonwound1 %>%
  PercentageFeatureSet(pattern = "^mt-")

obj_list$Nonwound1 %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3)

obj_list$Nonwound1 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")
obj_list$Nonwound1 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

obj_list$Nonwound1 <-
  obj_list$Nonwound1 %>%
  subset(nFeature_RNA > 200 &
           nFeature_RNA < 2500 &
           nCount_RNA < 8000 &
           percent.mt < 1)

obj_list$Nonwound1$sample <- "nonwound1"

# nonwound 2 qc -----------------------------------------------------------

obj_list$Nonwound2[["percent.mt"]] <- obj_list$Nonwound2 %>%
  PercentageFeatureSet(pattern = "^mt-")

obj_list$Nonwound2 %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3)

obj_list$Nonwound2 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")
obj_list$Nonwound2 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

obj_list$Nonwound2 <-
  obj_list$Nonwound2 %>%
  subset(nFeature_RNA > 200 &
           nFeature_RNA < 2200 &
           nCount_RNA < 6000 &
           percent.mt < 1)

obj_list$Nonwound2$sample <- "nonwound2"

# wound 1 qc -----------------------------------------------------------

obj_list$Wound1[["percent.mt"]] <- obj_list$Wound1 %>%
  PercentageFeatureSet(pattern = "^mt-")

obj_list$Wound1 %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3)

obj_list$Wound1 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")
obj_list$Wound1 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

obj_list$Wound1 <-
  obj_list$Wound1 %>%
  subset(nFeature_RNA > 200 &
           nFeature_RNA < 3000 &
           nCount_RNA < 6000 &
           percent.mt < 1.5)

obj_list$Wound1$sample <- "wound1"

# wound 2 qc -----------------------------------------------------------

obj_list$Wound2[["percent.mt"]] <- obj_list$Wound2 %>%
  PercentageFeatureSet(pattern = "^mt-")

obj_list$Wound2 %>%
  VlnPlot(features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3)

obj_list$Wound2 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "percent.mt")
obj_list$Wound2 %>%
  FeatureScatter(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

obj_list$Wound2 <-
  obj_list$Wound2 %>%
  subset(nFeature_RNA > 200 &
           nFeature_RNA < 1750 &
           nCount_RNA < 2250 &
           percent.mt < 1.5)

obj_list$Wound2$sample <- "wound2"

# SCTransform -------------------------------------------------------------
obj_list <- obj_list %>%
  lapply(
    FUN = function(x) {
      SCTransform(x, vst.flavor = "v2")
    }
  )

# Integration -------------------------------------------------------------
features <- obj_list %>%
  SelectIntegrationFeatures(nfeatures = 3000)
obj_list <- obj_list %>%
  PrepSCTIntegration(anchor.features = features)
anchors <- obj_list %>% # ~20 min with m1
  FindIntegrationAnchors(normalization.method = "SCT",
                         anchor.features = features)

combined_sct <- anchors %>%
  IntegrateData(normalization.method = "SCT")

combined_sct <- combined_sct %>%
  RunPCA(npcs = 30) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) |>
  FindClusters()

new.cluster.ids <- c("Fibroblasts 1", "Fibroblasts 2", "Fibroblasts 3", "Keratinocytes", "Fibroblasts 4", "Monocytes",
                     "Fibroblasts 5", "Macrophages", "N/A","Endothelial","Errector Pilli","Endothelial","Macrophages",
                     "Fibroblasts 6","Sk Mucle 1","SK Muscle 2","Adipocyte","Sk Muscle 3","Fibroblasts 7","Unknown 1", "unknown 2")
names(new.cluster.ids) <- levels(combined_sct)
combined_sct <- RenameIdents(combined_sct, new.cluster.ids)
DimPlot(combined_sct, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


combined_sct %>%
  write_rds("~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-integrated.rds")

combined_sct %>%
  DimPlot(label = T,
          label.box = T,
          repel = T) +
  theme(legend.position = "none")


ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/sctransform-integration-comparison/integrated/clusters-umap.png"
)

combined_sct %>%
  DimPlot(group.by = "sample")

ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/sctransform-integration-comparison/integrated/sample-umap.png"
)

combined_sct$type <-
  combined_sct$sample %>%
  gsub('[[:digit:]]+', '', .)

combined_sct %>%
  DimPlot(group.by = "type")

ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/sctransform-integration-comparison/integrated/injury-umap.png"
)
