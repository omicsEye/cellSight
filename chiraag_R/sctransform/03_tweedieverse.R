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
  readRDS("~/Box/snRNA_Cell_project/objects/2022-08-23_integrated-clustered.rds")

seur_obj$Celltype <- Idents(seur_obj)
seur_obj$type <- seur_obj$type %>%
  as.factor()
###For chekcing ####
i <- "Immune cells"
####################
for (i in seur_obj$Celltype %>% unique()) {
  print(i)
  obj_sub <- seur_obj %>%
    subset(Celltype == i,
           features = VariableFeatures(.))
  
  counts <- obj_sub %>%
    GetAssayData(assay = "RNA", slot = "counts")
  
  
  
  
  
  input_features <- counts |>
    as.matrix() |>
    t() |>
    as.data.frame()
  
  test <- Tweedieverse(
    input_features,
    obj_sub@meta.data[6] ,
    output = paste0('~/Desktop/Single_Cell_output/tweedieverse/cluster_', i),
    entropy_threshold = 0.0,
    base_model = 'CPLM',
    plot_heatmap = T,
    plot_scatter = T,
    reference = c("type, Normal"),
    cores = 7
  )
}

