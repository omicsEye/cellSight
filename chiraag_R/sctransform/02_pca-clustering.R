# load libraries ----------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clustifyr")

install.packages("clustree")
library(Seurat)
library(ggplot2)
library(clustree)
library(tidyverse)
library(clustifyr)
library(clustifyrdata)


## for MAC ###
seur_obj <-
  readRDS("~/Box/snRNA_Cell_project/objects/integrated.rds")

## For Windows ###

seur_obj <-
  readRDS("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/objects/sc-integrated.rds")



seur_obj <- seur_obj |>
  RunPCA()

DimPlot(seur_obj, group.by = "type")

DimPlot(seur_obj, group.by = "sample")

ElbowPlot(seur_obj, ndims = 50)

DefaultAssay(seur_obj) <- "integrated"
seur_obj <- seur_obj |>
  RunUMAP(dims = 1:11) |>
  FindNeighbors(dims = 1:11) |>
  FindClusters(resolution = c(.2, .4, .6, .8))

DimPlot(seur_obj,
        group.by = "integrated_snn_res.0.2",
        label = T,
        label.box = T) +
  theme(legend.position = "none")

DimPlot(seur_obj,
        group.by = "integrated_snn_res.0.4",
        label = T,
        label.box = T) +
  theme(legend.position = "none")

DimPlot(seur_obj,
        group.by = "integrated_snn_res.0.6",
        label = T,
        label.box = T) +
  theme(legend.position = "none")

seur_obj |>
  DimPlot(split.by = "sample")

Idents(seur_obj) <- "integrated_snn_res.0.6"
seur_obj <- seur_obj |>
  NormalizeData(assay = "RNA")

seur_obj$type <-
  seur_obj$sample %>%
  gsub('[[:digit:]]+', '', .)

seur_obj |>
  saveRDS(
    "~/Box/snRNA_Cell_project/objects/2022-08-23_integrated-clustered.rds"
    )

for (i in seur_obj$integrated_snn_res.0.6 |> unique()) {
  seur_obj |>
    FindConservedMarkers(ident.1 = i,
                         grouping.var = "type") |>
    write.csv(
      paste0(
        "~/Box/snRNA_Cell_project/analysis/chiraag-analysis/conserved-",
        i,
        ".csv"
      )
    )
}


