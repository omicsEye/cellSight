library(Seurat)
library(readr)

data <-
  readRDS("~/Box/snRNA_Cell_project/objects/2022-08-23_integrated-clustered.rds")

table(data$integrated_snn_res.0.6, data$sample) |>
  prop.table(margin = 1) |>
  write.csv(
    "~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/sample-cluster-metrics.csv"
  )


data <- data %>%
  NormalizeData(assay = "RNA")

Idents(data) <- "integrated_snn_res.0.6"

data |>
  FindAllMarkers(logfc.threshold = 0, 
                 assay = "RNA", 
                 only.pos = T) |>
  write.csv("~/Box/snRNA_CellRanger_Wound_nonWound/markers/sctransform-diff_qc-integration/top-markers.csv")
