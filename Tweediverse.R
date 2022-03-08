install.packages("devtools")
library(devtools)
install_github('bioBakery/Maaslin2', force = TRUE)
install_github("himelmallick/Tweedieverse", force = TRUE)
install_github('omicsEye/omicsArt', force = TRUE)
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
library(readr)

#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound")
#import the data from box
setwd("/Users/Rano/Desktop/Single_cell_output/")

####Read omics data table ##########################
data <- read.delim(
  'master_data.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
data <- test
#Adding the sample type

#read metadata
metadata <- read.delim(
  'master_metadata.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)
cell <- "Fibroblasts"
for (cell in unique(metadata$Cell_type)){
  cell_metadata <- metadata[which (metadata[,"Cell_type"] == cell),]
  cell_metadata <- cell_metadata[,c('Cell_type', 'Sample_type')]
  cell_metadata$Cell_type<-as.factor(cell_metadata$Cell_type)
  cell_metadata$Cell_type<-factor(cell_metadata$Cell_type, levels = c(cell))
  interim <- rownames(cell_metadata)
  cell_data <- subset(data, rownames(data)%in% interim)
  cell_data <- cell_data[, colSums(cell_data != 0) > 0]
  cell_metadata <- cell_metadata[,c( 'Sample_type')]
  #cell_metadata$['Sample_type']<-as.factor(cell_metadata$['Sample_type'])
  #cell_metadata$Sample_type<-factor(cell_metadata$['Sample_type'], levels = c("Wound","Nonwound"))
  #cell_data[is.na(cell_data)] <- 0
  if (TRUE){
    print(cell_metadata)
    Tweedieverse::Tweedieverse(cell_data,
                               cell_metadata,
                               paste0('analysis/Tweedieverse_output_', cell, sep = ""),
                               max_significance = 0.05,
                               base_model = "CPLM",
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               reference = NULL
          )
  }
}

ordplots(data = data, metadata = metadata, output = "/analysis/PCoA.pdf")
fakepcl <- list(meta=metadata, x=as.matrix(data),
                ns=dim(data)[2], nf=dim(data)[1])
heat_plot <- omicsArt:::pcl.heatmap(fakepcl, sqrtspace = T, gamma = 4, meta= T, show_colnames = F, show_rownames = T)

