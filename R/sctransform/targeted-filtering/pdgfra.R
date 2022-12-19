# load libraries ----------------------------------------------------------

library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(stringr)
library(Tweedieverse)

box_dir <- "~/Box/snRNA_CellRanger_Wound_nonWound/data/"
dirs <- list.dirs(path = box_dir,
                  recursive = F,
                  full.names = F)

obj_list <- dirs %>%
  set_names() %>%
  map(
    .f = function(x) {
      data <- paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
        Read10X() %>%
        CreateSeuratObject(min.cells = 3, min.features = 200) %>%
        AddMetaData(PercentageFeatureSet(., pattern = "^mt-"),
                    "percent.mt") %>%
        subset(nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt < 5)
      
      adgre1_expression <-
        GetAssayData(object = data,
                     assay = "RNA",
                     slot = "counts")["Pdgfra", ]
      pos_ids <- names(which(adgre1_expression > 0))
      data <- data %>%
        subset(cells = pos_ids)
      data$sample <- x
      data %>%
        SCTransform(vst.flavor = "v2")
    }
  )

features <- obj_list %>%
  SelectIntegrationFeatures(nfeatures = 3000)
obj_list <- obj_list %>%
  PrepSCTIntegration(anchor.features = features)

anchors <- obj_list %>%
  FindIntegrationAnchors(normalization.method = "SCT",
                         anchor.features = features)
combined_sct <- anchors %>%
  IntegrateData(normalization.method = "SCT", k.weight = 46)

combined_sct$injury <- str_remove(combined_sct$sample, "\\d")

combined_sct <- combined_sct %>%
  RunPCA(npcs = 100) %>%
  RunUMAP(dims = 1:50) %>%
  FindNeighbors(dims = 1:50) %>%
  FindClusters()

DimPlot(combined_sct)
ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/targeted filtering/sctransform/pdgfra/umap.png"
)
DimPlot(combined_sct, split.by = "injury")
ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/targeted filtering/sctransform/pdgfra/umap-split.png"
)

combined_sct <- combined_sct %>%
  NormalizeData(assay = "RNA")

combined_sct %>%
  FindAllMarkers(assay = "RNA",
                 only.pos = T) %>%
  write.csv(
    "~/Box/snRNA_CellRanger_Wound_nonWound/targeted filtering/sctransform/pdgfra/top-markers.csv"
  )


# tweedieverse ------------------------------------------------------------

DefaultAssay(combined_sct) <- "RNA"

for (i in combined_sct$seurat_clusters %>% unique()) {
  combined_sub <- combined_sct %>%
    subset(seurat_clusters == i)
  
  if (combined_sub$injury %>% table() %>% length() == 1) {
    next
  }

  counts <- combined_sub %>%
    GetAssayData(assay = "RNA", slot = "counts") %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
  
  Tweedieverse(
    counts,
    combined_sub@meta.data,
    output = paste0(
      '~/Box/snRNA_CellRanger_Wound_nonWound/targeted filtering/sctransform/pdgfra/tweedieverse/cluster-',
      i
    ),
    # Assuming demo_output exists
    fixed_effects = c('injury'),
    base_model = 'CPLM',
    adjust_offset = TRUE
  ) 
}
