library(Seurat)
library(tidyverse)
library(patchwork)

int <-
  read_rds("~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-integrated.rds")

nonint <-
  read_rds(
    "~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-non-integrated.rds"
  )

p1 <- int %>%
  DimPlot(group.by = "sample") +
  labs(title = "Integrated")
p2 <- nonint %>%
  DimPlot(group.by = "sample") +
  labs(title = "Non-integrated")

p1 + p2 +
  plot_layout(guides = "collect",)

ggsave(
  filename = "~/Box/snRNA_CellRanger_Wound_nonWound/sctransform-diff_qc-integration/sctransform-integration-comparison/split-comparision-plot.png",
  height = 6,
  width = 10,
  units = "in",
  dpi = "retina"
)

nw1 <- int %>%
  subset(sample == "nonwound1") %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  theme(legend.position = "none") +
  labs(title = "1D PNW1")

nw2 <- int %>%
  subset(sample == "nonwound2") %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  theme(legend.position = "none") +
  labs(title = "1D PNW2")

w1 <- int %>%
  subset(sample == "wound1") %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  theme(legend.position = "none") +
  labs(title = "1D PW1")

w2 <- int %>%
  subset(sample == "wound2") %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  theme(legend.position = "none") +
  labs(title = "1D PW2")

DefaultAssay(int) <- "SCT"


(nw1 + nw2) /
  (w1 + w2) +
  
  
  nonint$type <-
  nonint$sample %>%
  gsub('[[:digit:]]+', '', .)

int$type <-
  int$sample %>%
  gsub('[[:digit:]]+', '', .)

ggsave(
  "~/Box/snRNA_CellRanger_Wound_nonWound/analysis/sctransform/split-plot.png",
  dpi = "retina"
)

nonint_sample <- nonint %>%
  DimPlot(group.by = "sample") +
  labs(title = "Before integration") +
  theme(axis.text = element_blank())

nonint_injury <- nonint %>%
  DimPlot(group.by = "type") +
  labs(title = "Before integration") +
  theme(axis.text = element_blank())

int_sample <- int %>%
  DimPlot(group.by = "sample") +
  labs(title = "After integration") +
  theme(axis.text = element_blank())

int_injury <- int %>%
  DimPlot(group.by = "type") +
  labs(title = "After integration") +
  theme(axis.text = element_blank())

sc <- int %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  labs(title = NULL, x = NULL, y = NULL) + 
  theme(legend.position = "none") +
  theme(axis.text = element_blank())


x <- (nonint_sample + int_sample) +
  plot_layout(guides = "collect")

y <- nonint_injury + int_injury +
  plot_layout(guides = "collect")

x / y / sc +
  plot_layout(heights = c(1, 1, 2))

int %>%
  RidgePlot(features = c("Pdgfrb"),
            group.by = "sample")


Idents(int) <- "sample"

int <- int %>%
  NormalizeData(assay = "RNA")
