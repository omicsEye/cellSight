# or 
BiocManager::install("flowCore")
devtools::install_github("RyanYip-Kat/yipCat")
install.packages("yipCat_1.0.3.tar.gz")
# or 
devtools::install_github("RyanYip-Kat/yipCat",dependencies = TRUE)
BiocManager::install("ComplexHeatmap")
BiocManager::install("SummarizedExperiment")
BiocManager::install("DropletUtils")
BiocManager::install("Rsamtools")
library(yipCat)
setwd("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/objects")
seur_obj <-readRDS("scintegrated_final.rds") 
wounded <- subset(x = seur_obj, subset = (type  == 'Injured'))
naive_full <- subset(x = seur_obj, subset = (type  == 'Naive'))
cell_barcode = colnames(naive@assays$RNA@data)
features = rownames(naive@assays$RNA@data)
exportCellPhoneDB(naive,cells=cell_barcode,features=features,selectCol="Celltype",path2CPDB = '/Users/ranojoychatterjee/miniconda3/envs/cpdb',runCPDB=TRUE)  # selectCol : which column to caculate cell-cell interaction,runCPDB=TRUE,run cellphonedb backgroup
CPDBDotplot  #  cellphonedb result dotplot
CPDBHeatmaps  # cellphonedb result heatmap
remotes::install_github('lima1/sttkit')
library(sttkit)
cellphone_for_seurat(naive,org.Mm.eg.db,prefix = "results" )
import_cellphone(naive,cellphone_outpath = "results",org.Mm.eg.db)
matrix1 <- as.data.frame(as.matrix(seur_obj@assays$RNA@data))
matrix1 <- matrix1[rowSums(matrix1[,2:dim(matrix1)[2]])!=0,]

count_raw_meta <- GetAssayData(object = naive, slot = "counts")[,colnames(x = naive)]


human = useEnsembl(biomart='ensembl', dataset="hsapiens_gene_ensembl", mirror = "useast",version = 105)
mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",mirror = "useast",version = 105)
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(naive_full@assays$RNA@data) , mart = mouse, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)
print(head(genesV2))
matrix1 <- matrix1[match(genesV2$MGI.symbol,rownames(seur_obj),nomatch=F),]
count_norm_meta <- apply(count_raw_meta, 2, function(x) (x/sum(x))*10000)
matrixA <- count_norm_meta[match(genesV2$MGI.symbol,rownames(count_norm_meta),nomatch=F),]
matrix1$gene <- genesV2$Gene.stable.ID
Naive <- grepl('Naive',seur_obj@meta.data$type)
Injured <- grepl('Injured',seur_obj@meta.data$type)
Naive[match('gene',colnames(matrix1))] <- TRUE
Injured[match('gene',colnames(matrix1))] <- TRUE


write.table(matrix1[,Naive], 'Naive_filtered_hcount.csv',row.names=T,sep=',')
write.table(matrix1[,Injured], 'Injured_filtered_hcount.csv',row.names=T,sep=',')


metadata <- data.frame(cells=rownames(seur_obj@meta.data[grepl('Naive',seur_obj@meta.data$type),]),cluster=seur_obj@meta.data$Celltype[grepl('Naive',seur_obj@meta.data$type)])
metadata_s2 <- data.frame(cells=rownames(seur_obj@meta.data[!grepl('Naive',seur_obj@meta.data$type),]),cluster=seur_obj@meta.data$Celltype[!grepl('Naive',seur_obj@meta.data$type)]) ## Just negate grepl('state1',alldata@meta.data$stim),]
print('Writing Metadata')
write.csv(metadata, 'Naive_filtered_meta.csv', row.names=FALSE)
write.csv(metadata_s2, 'Injured_filtered_meta.csv', row.names=FALSE)





naive = naive_full@assays$RNA@data
meta = cbind(cell = rownames(naive_full@meta.data), cell_type=naive_full@meta.data$Celltype)

goi = unique(genesV2[,4])
keep_idx2 = match(goi, genesV2[,4])
naive = naive[keep_idx2, ]
genesV2 = genesV2[keep_idx2, ]
rownames(naive) = genesV2[, "Gene.stable.ID"]
matched_meta = meta[match(colnames(naive),meta[,1]),]

write.table(naive, "Naive_count.txt", sep ="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(matched_meta, "Naive_meta.txt", sep ="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

injured_match = wounded@assays$RNA@data
meta_inj = cbind(cell = rownames(wounded@meta.data), cell_type=wounded@meta.data$Celltype)

goi = unique(genesV2[,4])
keep_idx2 = match(goi, genesV2[,4])
injured_match = injured_match[keep_idx2, ]
genesV2 = genesV2[keep_idx2, ]
rownames(injured_match) = genesV2[, "Gene.stable.ID"]
matched_meta_inj = meta_inj[match(colnames(injured_match),meta_inj[,1]),]

write.table(injured_match, "Wounded_count.txt", sep ="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(matched_meta_inj, "Wounded_meta.txt", sep ="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



df<- read.csv("Naive_count.txt",sep = "\t")
library(tibble)
df <- tibble::rownames_to_column(df, "Gene")
write.table(df_try, "Naive_count.txt", sep ="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
df_1 <- read.csv("test_counts.txt",sep = "\t")
df_2 <- read.csv("Naive_meta.txt",sep = "\t")



BiocManager::install("InterCellar")
library(InterCellar)
InterCellar::run_app( reproducible = TRUE )



