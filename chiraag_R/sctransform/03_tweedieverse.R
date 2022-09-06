library(Seurat)
library(Tweedieverse)
library(tidyverse)

seur_obj <-
  readRDS("~/Box/snRNA_Cell_project/objects/2022-08-23_integrated-clustered.rds")

seur_obj$type <- seur_obj$type %>%
  as.factor()

for (i in seur_obj$integrated_snn_res.0.6 %>% unique()) {
  obj_sub <- seur_obj %>%
    subset(integrated_snn_res.0.6 == i,
           features = VariableFeatures(.))
  
  counts <- obj_sub %>%
    GetAssayData(assay = "RNA",
                 slot = "counts")
  
  input_features <- counts |>
    as.matrix() |>
    t() |>
    as.data.frame()
  
  test <- Tweedieverse(
    input_features,
    obj_sub@meta.data,
    output = paste0('chiraag_R/tweedieverse_output/cluster_', i),
    fixed_effects = c('type'),
    base_model = 'CPLM',
    adjust_offset = TRUE,
    cores = 4
  )
}
