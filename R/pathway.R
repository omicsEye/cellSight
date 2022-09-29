library(devtools)
#devtools::install_github('omicsEye/omepath', force = TRUE)
devtools::install_github('omicsEye/omePath', force = TRUE)
library(omePath)
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(KEGGREST)
library(tidyverse)
install.packages("msigdbr")
library(msigdbr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
library(fgsea)


setwd("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/")
data_dir <- "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/analysis/Tweedieverse_output_Pdgfrb/"
list.files(data_dir)
df <- read.table(file = 'analysis/Tweedieverse_output_Pdgfrb/all_results.tsv', sep = '\t', header = TRUE)
df$feature <- (df$feature)
#####Using Mgsidb######
all_gene_sets = msigdbr(species = "Mus musculus")
head(all_gene_sets)

original_gene_list <- df$coef
names(original_gene_list) <- (df$feature)

gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#This is the databse for mouse but it is only mapping to 5 genes
organism = "org.Mm.eg.db"

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

mapper_file<-gse@result

mapper_file<-mapper_file %>% mutate(core_enrichment = strsplit(as.character(core_enrichment), "/")) %>% unnest(core_enrichment)

names(mapper_file)[12] <- 'feature'

x <- enrichPathway(gene=dedup_ids$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
head(x)
dedup_ids = sort(dedup_ids$ENTREZID, decreasing = TRUE)
data(geneList, package="DOSE"
)

y <- gsePathway(kegg_gene_list, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = dedup_ids$ENTREZID)

organism = "org.Hs.eg.db"
#This is the databse for mouse but it is only mapping to 5 genes
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

original_gene_list <- df$coef
names(original_gene_list) <- toupper((df$feature))
df$feature <- toupper((df$feature))
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

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
kegg_organism = "mmu"
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               eps          = 0)
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

#emapplot(kk2)
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
#ridgeplot(gse) + labs(x = "enrichment distribution")

library(limma)
tab <- getGeneKEGGLinks(species="hsa")
tab$Symbol <- mapIds(org.Hs.eg.db, tab$GeneID,
                       column="SYMBOL", keytype="ENTREZID")
paths <- merge(df2,tab,by.x = "Y", by.y ="GeneID")

pathways <- getKEGGPathwayNames(species="hsa")
pathways = pathways[pathways$PathwayID %in% paths$PathwayID,]

input_data <- paths[,c("Y","PathwayID","coef","pval","qval")]
mapper_file <- pathways[,c("PathwayID","Description")]
colnames(mapper_file)[2] <- "Description"
input_data <- input_data %>% distinct(Y, .keep_all = TRUE)
rownames(input_data)<- input_data$Y
table1.df <- dplyr::inner_join(mapper_file,input_data, by="PathwayID")
setwd("C:/Users/ranoj/Box/Update_single_cell/")
omepath_results <- omePath(as.data.frame(df,row.names = df$feature),
                           "output_deepath",
                           as.data.frame(mapper_file,row.names = mapper_file$feature), 
                           pathway_col = "Description",
                           feature_col = "feature",
                           score_col = 'coef',
                           method = 'wilcox',
                           min_member = 2,
                           do_plot = TRUE)
