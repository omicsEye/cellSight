library(Seurat)
library(tidyverse)

data <-
  read_rds("~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-integrated.rds")

data$type <-
  data$sample %>%
  gsub('[[:digit:]]+', '', .)

data <- PrepSCTFindMarkers(data)

for (i in data$seurat_clusters %>% unique()) {
  data %>%
    FindConservedMarkers(ident.1 = i,
                         grouping.var = "type",
                         assay = "SCT") %>%
    write.csv(
      paste0(
        "~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/conserved-markers/cluster-",
        i,
        ".csv"
      )
    )
}

data %>%
  FindAllMarkers(assay = "SCT",
                 logfc.threshold = 0,
                 min.pct = 0,
                 only.pos = T) %>%
  write_csv("~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/top-markers.csv")

features <- c(
  "Adgre1",
  "Ptprc",
  "Krt10",
  "Plin1",
  "Acta2",
  "Itga8",
  "Col1a1",
  "Krt14",
  "Pdgfra",
  "Cd34"
)

DefaultAssay(data) <- "SCT"

features %>%
  walk(function(x) {
    data %>%
      FeaturePlot(x, max.cutoff = 3)
    ggsave(
      paste0(
        "~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/",
        x,
        "-umap.png"
      )
    )
  })

data %>%
  DimPlot(label = T,
          label.box = T,
          repel = T) +
  theme(legend.position = "none")

ggsave("~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/clusters-umap.png")

data %>%
  DimPlot(group.by = "sample")

ggsave("~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/sample-umap.png")

data %>%
  DimPlot(group.by = "type")

ggsave("~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/injury-umap.png")

Idents(data) <- "seurat_clusters"

data <- data %>%
  NormalizeData(assay = "RNA")

for (i in data$seurat_clusters %>% unique()) {
  data %>%
    FindMarkers(ident.1 = "wound", 
                group.by = "type", 
                subset.ident = i, 
                recorrect_umi = FALSE,
                logfc.threshold = 0,
                min.pct = 0,
                assay = "RNA") %>%
    write.csv(
      paste0(
        "~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/de/cluster-",
        i,
        ".csv"
      )
    )
}


cc_injury <- c("Ccl2", "Ccl3", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9", "Ccl11",  "Ccl19", "Ccl20", "Ccl25", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl4", "Cxcl5", "Cxcl12", "Cxcl13", "Cxcl14", "Csf1", "Csf2", "Csf3", "Fst", "Il11", "Il19", "Il33", "Osm", "Pgf", "Spp1", "Ppbp", "Tnfsf11")