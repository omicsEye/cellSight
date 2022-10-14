# load libraries ----------------------------------------------------------

install.packages("rlang")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("glmGamPoi")
install.packages("dplyr")  

library(dplyr)
library(glmGamPoi)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

box_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/"
dirs <- list.dirs(path = box_dir, 
                  recursive = F, 
                  full.names = F)

local_dir <- "C:/Users/ranoj/Desktop/Single_Cell_output/"

obj_list <- dirs %>%
  set_names() %>%
  map(.f = function(x) {
    paste0(box_dir, x, "/outs/filtered_feature_bc_matrix") %>%
      Read10X() %>%
      CreateSeuratObject()
  })

## Just changed the length to dimension to get both the size##

dim(obj_list$Nonwound1)
dim(obj_list$Nonwound2)
dim(obj_list$Wound1)
dim(obj_list$Wound2)

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
  lapply(FUN = function(x) {
    SCTransform(x, vst.flavor = "v2")
  })

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

combined_sct %>%
  saveRDS("~/Box/snRNA_Cell_project/objects/integrated.rds")