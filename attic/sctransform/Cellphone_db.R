library(org.Mm.eg.db)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(reticulate)
library(Seurat)
library(SeuratData)
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
devtools::install_github('satijalab/seurat-data')
BiocManager::install("scater",force = T)
library(scater)
install.packages("devtools")
devtools::install_github("whtns/seuratTools")
library(seuratTools)
setwd("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

seur_obj <-
  readRDS("objects/scintegrated_final.rds")

setwd("/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/")

Idents(seur_obj) = seur_obj$Celltype
BiocManager::install("biomaRt",force =T)
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(seur_obj@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))


SaveH5Seurat(seur_obj, filename = "seur_obj.h5Seurat")
Convert("seur_obj.h5Seurat", dest = "h5ad")
writeMM(seur_obj@assays$RNA@data, file = 'combined_counts_mtx/matrix.mtx')
# save gene and cell names
write(x = rownames(seur_obj@assays$RNA@data), file = "combined_counts_mtx/features.tsv")
write(x = colnames(seur_obj@assays$RNA@data), file = "combined_counts_mtx/barcodes.tsv")



table(seur_obj@meta.data$Celltype)

seur_obj@meta.data$Cell = rownames(seur_obj@meta.data)
df = seur_obj@meta.data[, c('Cell', 'Celltype')]
write.table(df, file ='combined_meta.tsv', sep = '\t', quote = F, row.names = F)


DEGs <- FindAllMarkers(seur_obj, logfc.threshold = 0,
                       min.pct = 0,
                       only.pos = T)

DEGs = DEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
write.table(DEGs, file ='combined_DEGs.tsv', sep = '\t', quote = F, row.names = F)


devtools::install_github("zhanghao-njmu/SCP")
library(SCP)
data("pancreas1k")
adata <- srt_to_adata(pancreas1k)
adata$write_h5ad("pancreas1k.h5ad")




seur_obj_sce <- as.SingleCellExperiment(seur_obj)



pbmc_ad <- Convert(from = seur_obj, to = "anndata", filename = "/convert.h5ad")
pbmc_ad



SaveH5Seurat(seur_obj, filename = "combined.h5Seurat")
Convert("combined.h5Seurat", dest = "h5ad")






###########Liana########################
remotes::install_github('saezlab/liana')
library(liana)
library(tidyverse)
library(OmnipathR)
library(liana)
library(purrr)
library(magrittr)

setwd("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/")

seur_obj <-
  readRDS("objects/scintegrated_final.rds")
DefaultAssay(seur_obj) <- "RNA"


rownames(seur_obj@assays$RNA@counts) <- stringr::str_to_title(rownames(seur_obj@assays$RNA@counts))
rownames(seur_obj@assays$RNA@data) <- stringr::str_to_title(rownames(seur_obj@assays$RNA@data))
show_homologene()
# Here, we will convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse

# Run LIANA with the orthologous resource
liana_res <- liana_wrap(seur_obj,
                        resource = 'custom', # resource has to be set to 'custom' to work with external resources
                        external_resource = ortholog_resource # provide orthologous resource

                        )

# aggregate
liana_res <- liana_res %>%
    liana_aggregate()

dplyr::glimpse(liana_res)
# Plot example
liana_res %>%
  liana_dotplot(ntop = 20)

write.csv(liana_res, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_combined_untrimmed.csv", row.names=FALSE, quote=FALSE) 

liana_trunc <- liana_res %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

p <- chord_freq(liana_trunc,cex =0.5) 

write.csv(liana_trunc, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_combined.csv", row.names=FALSE, quote=FALSE) 

############Wounded###############
DefaultAssay(wounded) <- "RNA"


rownames(wounded@assays$RNA@counts) <- stringr::str_to_title(rownames(wounded@assays$RNA@counts))
rownames(wounded@assays$RNA@data) <- stringr::str_to_title(rownames(wounded@assays$RNA@data))
show_homologene()
# Here, we will convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse

# Run LIANA with the orthologous resource
liana_res <- liana_wrap(wounded,
                        resource = 'custom', # resource has to be set to 'custom' to work with external resources
                        external_resource = ortholog_resource # provide orthologous resource
                        
)

# aggregate
liana_res <- liana_res %>%
  liana_aggregate()

dplyr::glimpse(liana_res)
# Plot example
liana_res %>%
  liana_dotplot(ntop = 20)

write.csv(liana_res, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_wounded_untrimmed.csv", row.names=FALSE, quote=FALSE) 

liana_trunc <- liana_res %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

p <- chord_freq(liana_trunc)
write.csv(liana_trunc, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_wounded.csv", row.names=FALSE, quote=FALSE) 

###########Naive##################
DefaultAssay(naive) <- "RNA"


rownames(naive@assays$RNA@counts) <- stringr::str_to_title(rownames(naive@assays$RNA@counts))
rownames(naive@assays$RNA@data) <- stringr::str_to_title(rownames(naive@assays$RNA@data))
show_homologene()
# Here, we will convert LIANA's Consensus resource to murine symbols
op_resource <- select_resource("Consensus")[[1]]

# Generate orthologous resource
ortholog_resource <- generate_homologs(op_resource = op_resource,
                                       target_organism = 10090) # mouse

liana_res <- liana_wrap(naive,
                        resource = 'custom', # resource has to be set to 'custom' to work with external resources
                        external_resource = ortholog_resource # provide orthologous resource
                        
)

# aggregate
liana_res <- liana_res %>%
  liana_aggregate()

dplyr::glimpse(liana_res)
# Plot example
liana_res %>%
  liana_dotplot(ntop = 20)

write.csv(liana_trunc, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_naive_untrimmed.csv", row.names=FALSE, quote=FALSE) 

liana_trunc <- liana_res %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

p <- chord_freq(liana_trunc)


p <- chord_freq(liana_trunc,cex =0.5) 

write.csv(liana_trunc, "/Users/ranojoychatterjee/Desktop/Single_cell_output/cellphonedb/ligand_receptor_naive.csv", row.names=FALSE, quote=FALSE) 

