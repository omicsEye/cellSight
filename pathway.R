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
df$feature <- toupper(df$feature)
organism = "org.Hs.eg.db"
#This is the databse for mouse but it is only mapping to 5 genes
#organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- df$coef
names(original_gene_list) <- toupper((df$feature))

ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$feature %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df$coef

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
#This is for mice
#kegg_organism = "mmu"
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

