library(Seurat)

seur_obj <- readRDS("objects/2022-08-23_integrated-clustered.rds")

Idents(seur_obj) <- "integrated_snn_res.0.6"

seur_obj <- seur_obj |>
  NormalizeData(assay = "RNA")

DefaultAssay(seur_obj) <- "RNA"

features <- c("Adgre1", "Ptprc","Krt10", "Plin1",
              "Acta2", "Itga8", "Col1a1", "Krt14",
              "Pdgfra", "Cd34")

seur_obj |>
  FeaturePlot(features)
