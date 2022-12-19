library(Seurat)
library(tidyverse)
library(patchwork)

int <-
  read_rds("~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-integrated.rds")

nonint <-
  read_rds(
    "~/Box/snRNA_CellRanger_Wound_nonWound/objects/2022-09-25_sc-non-integrated.rds"
  )

nonint$type <-
  nonint$sample %>%
  gsub('[[:digit:]]+', '', .)

int$type <-
  int$sample %>%
  gsub('[[:digit:]]+', '', .)


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

sc_int <- int %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.box = T,
    label.size = 3
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_blank())

sc_nonint <- non_int %>%
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

x / y / sc_int +
  plot_layout(heights = c(1, 1, 2))
