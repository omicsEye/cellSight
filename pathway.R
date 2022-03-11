library(devtools)
devtools::install_github('omicsEye/deepath', force = TRUE)
library(deepath)
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
data_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/analysis/Tweedieverse_output_Fibroblasts/"
list.files(data_dir)
df <- read.table(file = 'analysis/Tweedieverse_output_Fibroblasts/significant_results.tsv', sep = '\t', header = TRUE)
organism = "org.Hs.eg.db"
#This is the databse for mouse but it is only mapping to 5 genes
#organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- df$coef
names(original_gene_list) <- toupper((df$feature))

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
