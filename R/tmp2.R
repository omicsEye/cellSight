library(Seurat)
library(tidyverse)

seur_obj <- readRDS("~/Box/snRNA_Cell_project/objects/2022-08-23_integrated-clustered.rds")

Idents(seur_obj) <- "integrated_snn_res.0.6"

seur_obj <- seur_obj |>
  NormalizeData(assay = "RNA")

DefaultAssay(seur_obj) <- "RNA"

features <- c("Adgre1", "Ptprc","Krt10", "Plin1",
              "Acta2", "Itga8", "Col1a1", "Krt14",
              "Pdgfra", "Cd34")


for (feature in features) {
  seur_obj |>
    FeaturePlot(feature, split.by = "sample")
  ggsave(
    paste0(
      "~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/split-feature-plots/",
      feature,
      ".png"
    ))
}


seur_obj %>%
  subset(sample %in% c("wound1", "wound2")) %>%
  FeaturePlot("Adgre1", split.by = "sample", interactive = T)
  