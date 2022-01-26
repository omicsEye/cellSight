library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(janitor)


#setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
setwd("~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

#data_dir <- 'C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data'
data_dir <- '~/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/data'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

pbmc.data_Wound1 <- Read10X(data.dir = "data/Wound1/filtered_feature_bc_matrix")
pbmc_Wound1 <- CreateSeuratObject(counts = pbmc.data_Wound1, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc_Wound1

pbmc.data_Wound2 <- Read10X(data.dir = "data/Wound2/filtered_feature_bc_matrix")
pbmc_Wound2 <- CreateSeuratObject(counts = pbmc.data_Wound2, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc_Wound2

pbmc.data_Nonwound1 <- Read10X(data.dir = "data/Nonwound1/filtered_feature_bc_matrix")
pbmc_Nonwound1 <- CreateSeuratObject(counts = pbmc.data_Nonwound1, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc_Nonwound1

pbmc.data_Nonwound2 <- Read10X(data.dir = "data/Nonwound2/filtered_feature_bc_matrix")
pbmc_Nonwound2 <- CreateSeuratObject(counts = pbmc.data_Nonwound2, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc_Nonwound2

initial <- rownames(pbmc.markers)
wound1 <- as.data.frame(pbmc.data_Wound1)
wound2 <- as.data.frame(pbmc.data_Wound2)
unwound1 <- as.data.frame(pbmc.data_Nonwound1)
unwound2 <- as.data.frame(pbmc.data_Nonwound2)

#compare_df_rows_same(test, test_un)
#final <- rbind(test,test_un)
initial <- row.names(pbmc.markers)
wound1 <- subset(wound1, rownames(wound1) %in% initial)
wound2 <- subset(wound2, rownames(wound2) %in% initial)
unwound1 <- subset(unwound1, rownames(unwound1) %in% initial)
unwound2 <- subset(unwound2, rownames(unwound2) %in% initial)
final <- rownames(wound1)
rownames(wound1) <- paste(row.names(wound1), c("_s1"))
rownames(wound2) <- paste(row.names(wound2), c("_s2")) 
rownames(unwound1) <- paste(row.names(unwound1), c("_s3")) 
rownames(unwound2) <- paste(row.names(unwound2), c("_s4")) 
test <- dplyr::bind_rows(wound1, wound2)
test <- dplyr::bind_rows(test, unwound1)
col_name <- subset(pbmc.markers, rownames(pbmc.markers)%in% final)
pbmc.markers <- subset(pbmc.markers, rownames(pbmc.markers) %in% final)
name <- c(row_names)
df <- data.frame(name)
write.table(df,"meta_data.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(test,"data.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
