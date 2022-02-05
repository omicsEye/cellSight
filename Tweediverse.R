
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
library(stringr)
library(readxl)
library(dplyr)

setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound")


####Read omics data table ##########################
data <- read.delim(
  'analysis/data/data_Wound1.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

#read metadata
metadata <- read.delim(
  'analysis/data/meta_data_Wound1.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
#cell <- "B-cell"
for (cell in unique(metadata$Cell_type)){
  cell_metadata <- metadata[metadata$Cell_type == cell,]
  Tweedieverse::Tweedieverse(data,
                             cell_metadata,
                             paste('analysis/Tweedieverse_output_', cell),
                             max_significance = 0.05,
                             base_model = "CPLM",
                             plot_heatmap = T,
                             plot_scatter = T,
                             standardize = F,
                            reference = NULL
)
}

ordplots(data = data, metadata = metadata, output = "/analysis/PCoA.pdf")
fakepcl <- list(meta=metadata, x=as.matrix(data),
                ns=dim(data)[2], nf=dim(data)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 4, meta= T, show_colnames = F, show_rownames = T)

