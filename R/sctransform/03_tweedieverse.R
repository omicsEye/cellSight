install.packages("devtools")
install.packages(c('dplyr','pbapply', 'lme4', 'lmerTest', 'car', 'cplm', 'pscl', 'logging', 'ggrepel', 'gridExtra', 'future', 'cowplot', 'Hmisc', 'MultiDataSet', 'TSP', 'htmlTable', 'igraph', 'insight', 'lubridate', 'mgcv', 'mvtnorm', 'optparse', 'parameters', 'pillar', 'pkgload', 'plotly', 'rlang', 'rvest', 'seriation', 'usethis', 'viridis', 'signal', 'tsne', 'openxlsx', 'readxl', 'xfun', 'yulab.utils', "labdsv", "seriation","diffusionMap"), repos='http://cran.r-project.org')
BiocManager::install("ropls")
devtools::install_github('omicsEye/omicsArt', dependencies = TRUE,force = TRUE)
library(Seurat)
library(Tweedieverse)
library(tidyverse)
library(gdata)
library(pheatmap)
library(vegan)
library(corrplot)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(Hmisc)
library(Maaslin2)
library(Tweedieverse)
library(omicsArt)


seur_obj <-
  readRDS("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/object/final_integrated_clustered1118.rds")

seur_obj <-
  readRDS("C:/Users/ranoj/Desktop/Single_cell_output/objects/final_integrated_clustered1118.rds")


seur_obj <- combined_sct

seur_obj$Celltype <- Idents(seur_obj)
seur_obj$type <- seur_obj$type %>%
  as.factor()

saveRDS(seur_obj,"scintegrated_final.rds")

seur_obj <-
  readRDS("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/scintegrated_final.rds")
###For chekcing ####
i <- "Adipocyte"
####################
for (i in seur_obj$Celltype %>% unique()) {
  print(i)
  obj_sub <- seur_obj %>%
    subset(Celltype == i)
  
  counts <- obj_sub %>%
    GetAssayData(assay = "RNA", slot = "counts")
  
  
  
  
  
  input_features <- counts |>
    as.matrix() |>
    t() |>
    as.data.frame()
  
  test <- Tweedieverse(
    input_features,
    obj_sub@meta.data[6] ,
    output = paste0('C:/Users/ranoj/Desktop/Single_cell_output/tweedieverse/cluster_', i),
    prev_threshold = 0.0,
    entropy_threshold = 0.0,
    base_model = 'CPLM',
    plot_heatmap = T,
    plot_scatter = T,
    reference = c("type, Normal")
  )
}

wounded <- subset(x = seur_obj, subset = (type  == 'wound'))

#test.subset <- subset(x = epithelial, subset = (stim == "Healthy" | stim == "another_condition"))

