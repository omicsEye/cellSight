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
library(tidyverse)
install.packages("read.delim", dependencies=TRUE)
library(utils)

#setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound")
#import the data from box
setwd("C:/Users/ranoj/Desktop/Single_Cell_output/")

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
#data <- test
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

uni1 <- read.delim(
  'Nonwound_1_matched.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

uni2 <- read.delim(
  'Nonwound_2_matched.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

wou1 <- read.delim(
  'Wound_1_matched.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

wou2 <- read.delim(
  'Wound_2_matched.tsv',
  sep = '\t',
  header = TRUE,
  fill = FALSE,
  comment.char = "" ,
  check.names = FALSE,
  row.names = 1
)

rownames(uni1) <- gsub("Nonwound1_uninjuredskin_", "", rownames(uni1), fixed=TRUE)
rownames(uni1) <- paste(row.names(uni1), c("_s3"))
View(uni1)

rownames(uni2) <- gsub("Nonwound2_uninjuredskin_", "", rownames(uni2), fixed=TRUE)
rownames(uni2) <- paste(row.names(uni2), c("_s4"))

rownames(wou1) <- gsub("Wound1_injuredskin_", "", rownames(wou1), fixed=TRUE)
rownames(wou1) <- paste(row.names(wou1), c("_s1"))

rownames(wou2) <- gsub("Wound2_injuredskin_", "", rownames(wou2), fixed=TRUE)
rownames(wou2) <- paste(row.names(wou2), c("_s2"))

test <- dplyr::bind_rows(wou1, wou2)
test <- dplyr::bind_rows(test, uni1)
test <- dplyr::bind_rows(test, uni2)

final <- rownames(test)
data_matched <- subset(data, rownames(data)%in% final)
metadata_matched <- subset(metadata, rownames(metadata)%in% final)

data_unmatched <- subset(data, !(rownames(data)%in% final))
metadata_unmatched <- subset(metadata, !(rownames(metadata)%in% final))


cell <- "Stromal"
for (cell in unique(metadata_matched$Cell_type)){
  cell_metadata <- metadata_matched[(metadata_matched[,"Cell_type"] == cell),]
  cell_metadata <- cell_metadata[,c('Cell_type','Sample_type')]
  cell_metadata$Cell_type<-as.factor(cell_metadata$Cell_type)
  cell_metadata$Cell_type<-factor(cell_metadata$Cell_type, levels = c(cell))
  interim_meta <- rownames(cell_metadata)
  cell_data <- subset(data_matched, rownames(data_matched)%in% interim_meta)
  interim_data <- rownames(cell_data)
  cell_metadata <- subset(cell_metadata, rownames(cell_metadata)%in% interim_data)

  cell_data <- cell_data[, colSums(cell_data != 0) > 0]
  cell_metadata <- data.frame(cell_metadata[,c( 'Sample_type')])
  #cell_metadata <- data.frame(cell_metadata)
  names(cell_metadata)<-'Sample_type'
  rownames(cell_metadata) = rownames(cell_data)
  cell_metadata$'Sample_type'<- as.factor(cell_metadata$'Sample_type')
  cell_metadata$Sample_type<-factor(cell_metadata$'Sample_type', levels = c("Nonwound","Wound"))
  cell_data[is.na(cell_data)] <- 0
  
  # cell_metadata <- metadata[which (metadata[,"Cell_type"] == cell),]
  # cell_metadata <- cell_metadata[,c('Cell_type','Sample_type')]
  # cell_metadata$Cell_type<-as.factor(cell_metadata$Cell_type)
  # cell_metadata$Cell_type<-factor(cell_metadata$Cell_type, levels = c(cell))
  # interim <- rownames(cell_metadata)
  # cell_data <- subset(data, rownames(data)%in% interim)
  # cell_data <- cell_data[, colSums(cell_data != 0) > 0]
  # cell_metadata <- cell_metadata[,c( 'Sample_type')]
  # #cell_metadata$['Sample_type']<-as.factor(cell_metadata$['Sample_type'])
  # #cell_metadata$Sample_type<-factor(cell_metadata$['Sample_type'], levels = c("Wound","Nonwound"))
  # #cell_data[is.na(cell_data)] <- 0
  if (TRUE){
    print(cell_metadata)
    Tweedieverse::Tweedieverse(cell_data,
                               cell_metadata,
                               paste0('analysis/Tweedieverse_output_', cell, sep = ""),
                               entropy_threshold = 0.0,
                               base_model = 'CPLM',
                               plot_heatmap = T,
                               plot_scatter = T,
                               standardize = F,
                               reference = c("Sample_type,Nonwound")
          )
  }
}



cell <- "Stromal"
for (cell in unique(metadata_unmatched$Cell_type)){
  cell_metadata <- metadata_unmatched[which (metadata_unmatched[,"Cell_type"] == cell),]
  cell_metadata <- cell_metadata[,c('Cell_type','Sample_type')]
  cell_metadata$Cell_type<-as.factor(cell_metadata$Cell_type)
  cell_metadata$Cell_type<-factor(cell_metadata$Cell_type, levels = c(cell))
  interim <- rownames(cell_metadata)
  cell_data <- subset(data_unmatched, rownames(data_unmatched)%in% interim)
  cell_data <- cell_data[, colSums(cell_data != 0) > 0]
  #cell_metadata <- cell_metadata[,c( 'Sample_type')]
  #cell_metadata$['Sample_type']<-as.factor(cell_metadata$['Sample_type'])
  #cell_metadata$Sample_type<-factor(cell_metadata$['Sample_type'], levels = c("Wound","Nonwound"))
  #cell_data[is.na(cell_data)] <- 0
  
  cell_metadata <- data.frame(cell_metadata[,c( 'Sample_type')])
  #cell_metadata <- data.frame(cell_metadata)
  names(cell_metadata)<-'Sample_type'
  rownames(cell_metadata) = rownames(cell_data)
  cell_metadata$'Sample_type'<- as.factor(cell_metadata$'Sample_type')
  cell_metadata$Sample_type<-factor(cell_metadata$'Sample_type', levels = c("Wound","Nonwound"))
  cell_data[is.na(cell_data)] <- 0
  # cell_metadata <- metadata[which (metadata[,"Cell_type"] == cell),]
  # cell_metadata <- cell_metadata[,c('Cell_type','Sample_type')]
  # cell_metadata$Cell_type<-as.factor(cell_metadata$Cell_type)
  # cell_metadata$Cell_type<-factor(cell_metadata$Cell_type, levels = c(cell))
  # interim <- rownames(cell_metadata)
  # cell_data <- subset(data, rownames(data)%in% interim)
  # cell_data <- cell_data[, colSums(cell_data != 0) > 0]
  # cell_metadata <- cell_metadata[,c( 'Sample_type')]
  # #cell_metadata$['Sample_type']<-as.factor(cell_metadata$['Sample_type'])
  # #cell_metadata$Sample_type<-factor(cell_metadata$['Sample_type'], levels = c("Wound","Nonwound"))
  # #cell_data[is.na(cell_data)] <- 0
  if (TRUE){
    print(cell_metadata)
    Tweedieverse::Tweedieverse(cell_data,
                               cell_metadata,
                               paste0('analysis/unmatched/Tweedieverse_output_', cell, sep = ""),
                               entropy_threshold = 0.0,
                               base_model = 'CPLM',
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

Maaslin2::Maaslin2(
  input_data  = cell_data,
  input_metadata =  cell_metadata,
  output =  paste0('analysis/Maaslin_output_', cell, sep = ""),
  max_significance= 0.1,
  analysis_method = 'LM',
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = FALSE
)

