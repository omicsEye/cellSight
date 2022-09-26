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
  plot_layout(guides = "collect", ) 

ggsave(
  filename = "~/Box/snRNA_CellRanger_Wound_nonWound/sctransform-diff_qc-integration/sctransform-integration-comparison/split-comparision-plot.png",
  height = 6,
  width = 10,
  units = "in",
  dpi = "retina"
)
