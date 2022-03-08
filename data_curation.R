library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SingleR)
library(janitor)
library(tidyr)


setwd("/Users/Rano/Desktop/Single_cell_output/")

data_dir <- '/Users/Rano/Desktop/Single_cell_outptut/'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


#Import the data from the BOX folder, reconfigure the first column as the row name
# and also add the desired suffix to differentiaite between the samples to rowname

#initial <- rownames(pbmc.markers)
wound1 <- read.table(file = 'analysis/data/data_labelled_Wound1.tsv', sep = '\t', header = TRUE)
rownames(wound1) <- wound1[,1]
wound1[,1] <- NULL
rownames(wound1) <- paste(row.names(wound1), c("_s1"))


wound2 <- read.table(file = 'analysis/data/data_labelled_Wound2.tsv', sep = '\t', header = TRUE)
rownames(wound2) <- wound2[,1]
wound2[,1] <- NULL
rownames(wound2) <- paste(row.names(wound2), c("_s2"))

nonwound1 <- read.table(file = 'analysis/data/data_labelled_Nonwound1.tsv', sep = '\t', header = TRUE)
rownames(nonwound1) <- nonwound1[,1]
nonwound1[,1] <- NULL
rownames(nonwound1) <- paste(row.names(nonwound1), c("_s3"))


nonwound2 <- read.table(file = 'analysis/data/data_labelled_Nonwound2.tsv', sep = '\t', header = TRUE)
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

meta_wound1 <- read.table(file = 'analysis/data/Label_Wound1.tsv', sep = '\t', header = TRUE)
rownames(meta_wound1) <- meta_wound1[,1]
meta_wound1[,1] <- NULL
rownames(meta_wound1) <- paste(row.names(meta_wound1), c("_s1"))
meta_wound1$Sample_type <- "Wound"
meta_wound1$Sample_number <- "Sample_1"
colmn <- paste("col", 1:3)
meta_wound1<- meta_wound1 %>% separate(ref.data.fine, sep = " ", into = colmn, remove = FALSE)
meta_wound1 <- meta_wound1[ , -which(names(meta_wound1) %in% c("col 2","col 3"))]
names(meta_wound1)[names(meta_wound1) == 'col 1'] <- "Cell_type"

meta_wound2 <- read.table(file = 'analysis/data/Label_Wound2.tsv', sep = '\t', header = TRUE)
rownames(meta_wound2) <- meta_wound2[,1]
meta_wound2[,1] <- NULL
rownames(meta_wound2) <- paste(row.names(meta_wound2), c("_s2"))
meta_wound2$Sample_type <- "Wound"
meta_wound2$Sample_number <- "Sample_2"
colmn <- paste("col", 1:3)
meta_wound2<- meta_wound2 %>% separate(ref.data.fine, sep = " ", into = colmn, remove = FALSE)
meta_wound2 <- meta_wound2[ , -which(names(meta_wound2) %in% c("col 2","col 3"))]
names(meta_wound2)[names(meta_wound2) == 'col 1'] <- "Cell_type"

meta_nonwound1 <- read.table(file = 'analysis/data/Label_Nonwound1.tsv', sep = '\t', header = TRUE)
rownames(meta_nonwound1) <- meta_nonwound1[,1]
meta_nonwound1[,1] <- NULL
rownames(meta_nonwound1) <- paste(row.names(meta_nonwound1), c("_s3"))
meta_nonwound1$Sample_type <- "Nonwound"
meta_nonwound1$Sample_number <- "Sample_3"
colmn <- paste("col", 1:3)
meta_nonwound1<- meta_nonwound1 %>% separate(ref.data.fine, sep = " ", into = colmn, remove = FALSE)
meta_nonwound1 <- meta_nonwound1[ , -which(names(meta_nonwound1) %in% c("col 2","col 3"))]
names(meta_nonwound1)[names(meta_nonwound1) == 'col 1'] <- "Cell_type"


meta_nonwound2 <- read.table(file = 'analysis/data/Label_Nonwound2.tsv', sep = '\t', header = TRUE)
rownames(meta_nonwound2) <- meta_nonwound2[,1]
meta_nonwound2[,1] <- NULL
rownames(meta_nonwound2) <- paste(row.names(meta_nonwound2), c("_s4"))
meta_nonwound2$Sample_type <- "Nonwound"
meta_nonwound2$Sample_number <- "Sample_4"
colmn <- paste("col", 1:3)
meta_nonwound2<- meta_nonwound2 %>% separate(ref.data.fine, sep = " ", into = colmn, remove = FALSE)
meta_nonwound2 <- meta_nonwound2[ , -which(names(meta_nonwound2) %in% c("col 2","col 3"))]
names(meta_nonwound2)[names(meta_nonwound2) == 'col 1'] <- "Cell_type"

##Merging the 4 seperate samples as a master meta dataset
meta_test <- dplyr::bind_rows(meta_wound1, meta_wound2)
meta_test <- dplyr::bind_rows(meta_test, meta_nonwound1)
meta_test <- dplyr::bind_rows(meta_test, meta_nonwound2)

##Removing the cell barcodes from the datas  that are not present in the metadata set
final <- rownames(meta_test)
final_data <- subset(test, rownames(test)%in% final)

write.table(meta_test,"master_metadata.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
write.table(test,"master_data.tsv",  sep = "\t", eol = "\n", quote = F, col.names = NA, row.names = T)
