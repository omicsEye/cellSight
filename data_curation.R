library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(janitor)


setwd("/Users/Rano/Desktop/Single_Cell_Wound/")

data_dir <- '/Users/Rano/Desktop/Single_Cell_Wound/'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


#Import the data from the BOX folder, reconfigure the first column as the row name
# and also add the desired suffix to differentiaite between the samples to rowname

#initial <- rownames(pbmc.markers)
wound1 <- read.table(file = 'analysis/data/data_Wound1.tsv', sep = '\t', header = TRUE)
rownames(wound1) <- wound1[,1]
wound1[,1] <- NULL
rownames(wound1) <- paste(row.names(wound1), c("_s1"))

wound2 <- read.table(file = 'analysis/data/data_Wound2.tsv', sep = '\t', header = TRUE)
rownames(wound2) <- wound2[,1]
wound2[,1] <- NULL
rownames(wound2) <- paste(row.names(wound2), c("_s2"))

nonwound1 <- read.table(file = 'analysis/data/data_Nonwound1.tsv', sep = '\t', header = TRUE)
rownames(nonwound1) <- nonwound1[,1]
nonwound1[,1] <- NULL
rownames(nonwound1) <- paste(row.names(nonwound1), c("_s3"))


nonwound2 <- read.table(file = 'analysis/data/data_Nonwound2.tsv', sep = '\t', header = TRUE)
rownames(nonwound2) <- nonwound2[,1]
nonwound2[,1] <- NULL
rownames(nonwound2) <- paste(row.names(nonwound2), c("_s4"))

#compare_df_rows_same(test, test_un)
#final <- rbind(test,test_un)
##Merging the 4 seperate samples as a master dataset
test <- dplyr::bind_rows(wound1, wound2)
test <- dplyr::bind_rows(test, nonwound1)
test <- dplyr::bind_rows(test, nonwound2)

##Import the metadata from the BOX folder

meta_wound1 <- read.table(file = 'analysis/data/meta_data_Wound1.tsv', sep = '\t', header = TRUE)
rownames(meta_wound1) <- meta_wound1[,1]
meta_wound1[,1] <- NULL
rownames(meta_wound1) <- paste(row.names(meta_wound1), c("_s1"))

meta_wound2 <- read.table(file = 'analysis/data/meta_data_Wound2.tsv', sep = '\t', header = TRUE)
rownames(meta_wound2) <- meta_wound2[,1]
meta_wound2[,1] <- NULL
rownames(meta_wound2) <- paste(row.names(meta_wound2), c("_s2"))

meta_nonwound1 <- read.table(file = 'analysis/data/meta_data_Nonwound1.tsv', sep = '\t', header = TRUE)
rownames(meta_nonwound1) <- meta_nonwound1[,1]
meta_nonwound1[,1] <- NULL
rownames(meta_nonwound1) <- paste(row.names(meta_nonwound1), c("_s3"))


meta_nonwound2 <- read.table(file = 'analysis/data/meta_data_Nonwound2.tsv', sep = '\t', header = TRUE)
rownames(meta_nonwound2) <- meta_nonwound2[,1]
meta_nonwound2[,1] <- NULL
rownames(meta_nonwound2) <- paste(row.names(meta_nonwound2), c("_s4"))

##Merging the 4 seperate samples as a master meta dataset
meta_test <- dplyr::bind_rows(meta_wound1, meta_wound2)
meta_test <- dplyr::bind_rows(meta_test, meta_nonwound1)
meta_test <- dplyr::bind_rows(meta_test, meta_nonwound2)

##Removing the cell barcodes from the datas  that are not present in the metadata set
final <- rownames(meta_test)
try <- subset(test, rownames(test)%in% final)

write.table(try,"master_metadata.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(meta_test,"master_data.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
