# libraries ---------------------------------------------------------------

library(devtools)
library(omePath)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(limma)
library(KEGGREST)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(org.Mm.eg.db)

data_dir <-
  setwd("/Users/ranojoychatterjee/Desktop/Single_cell_output/tweedieverse_rerun/cluster_Fibroblasts 7/")
df <- read.table(
  file = paste0(data_dir, "/all_results.tsv"),
  sep = '\t',
  header = TRUE
)

ids <-
  bitr(
    df$feature,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

tab <- getGeneKEGGLinks(species.KEGG = "mmu")
tab$Symbol <- mapIds(org.Mm.eg.db,
                     tab$GeneID,
                     column = "SYMBOL",
                     keytype = "ENTREZID")

gene_intersect <- intersect(tab$Symbol, df$feature)
df <- df %>%
  filter(feature %in% gene_intersect)

tab <- tab %>%
  filter(Symbol %in% gene_intersect)

omepath_results <- omePath(
  as.data.frame(df, row.names = df$feature),
  "output_deepath",
  as.data.frame(tab, row.names = tab$Symbol),
  pathway_col = "PathwayID",
  feature_col = "Symbol",
  score_col = 'coef',
  method = 'wilcox',
  min_member = 2,
  do_plot = TRUE
)

