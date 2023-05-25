library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
install.packages("svglite")
library(svglite)
my_cols = brewer.pal(4,"Dark2")
int <-
  read_rds("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/objects/sc-integrated.rds")
int <-
  read_rds("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/objects/sc-integrated.rds")
nonint <-
  read_rds(
    "C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/objects/sc-non-integrated.rds"
  )
nonint <-
  read_rds("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/objects/sc-non-integrated.rds")

seur_obj <- read_rds("/Users/ranojoychatterjee/Library/CloudStorage/Box-Box/snRNA_CellRanger_Wound_nonWound/objects/scintegrated_final.rds")

#int<- read_rds("scintegrated_final.rds")


nonint$type <-
  nonint$sample %>%
  gsub('[[:digit:]]+', '', .)

int$type <-
  int$sample %>%
  gsub('[[:digit:]]+', '', .)


nonint_sample <- nonint %>%
  DimPlot(group.by = "sample",cols=alpha(my_cols),pt.size=0.7,raster =T) +
  labs(title = "Before integration") +
  theme(axis.text = element_blank()) +
   NoAxes()+ggtitle(" ")+NoLegend()+
  scale_fill_brewer(palette = "Dark2")

nonint_injury <- nonint %>%
  DimPlot(group.by = "type",cols=alpha(my_cols),pt.size=0.7,raster =T) +
  labs(title = "Before integration") +
  theme(axis.text = element_blank()) +
  NoAxes()+ggtitle(" ")+NoLegend()+
  scale_fill_brewer(palette = "Dark2")

int_sample <- int %>%
  DimPlot(group.by = "sample",cols=alpha(my_cols),pt.size=0.7,raster =T) +
  labs(title = "After integration") +
  theme(axis.text = element_blank(),text = element_text(size = 8)) +
  NoAxes()+ggtitle(" ")+guides(color = guide_legend(override.aes = list(size=2), ncol=1) )+
  scale_fill_brewer(palette = "Dark2")

int_injury <- int %>%
  DimPlot(group.by = "type",raster = T) +
  labs(title = "After integration",cols=alpha(my_cols),pt.size=0.7) +
  theme(axis.text = element_blank(),text = element_text(size = 8)) +
  NoAxes()+ggtitle(" ")+ guides(color = guide_legend(override.aes = list(size=2), ncol=1) )+
  scale_fill_brewer(palette = "Dark2")

sc_int <- seur_obj %>%
  DimPlot(
    group.by = "Celltype",
    label = T,
    label.size = 2,
    raster = T, 
    pt.size = .01,
    repel =T,
    cols = c('KC' = "#49B500", 'SG 1' = '#00B4F0', 'SG 2' = '#FF61CC', 'Mono' = '#E36EF6',
                      'Mac' = '#FF64B0','Macro' ='#9F8CFF', 'Mast' = '#619CFF', 'FB 1' = '#CD9600',
                      'FB 2' = '#00BBDC','FB 3' = '#00C08B', 'FB 4' = '#00A9FF','FB 5' = '#9DA700','FB 6' = '#DE8C00','FB 7' = '#C77CFF',
                      'Adipo' = '#B79F00','Mu 1' = '#00BE67','Mu 2' = '#00BFC4', 'Mu 3' = '#00C1A9', 'APM' = '#00BA38', 'EC 1' = "#7CAE00",
                      ' EC 2' =  "#F564E3")
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_blank(),text = element_text(size = 2)) +
  NoAxes()
  
ggsave("/Users/ranojoychatterjee/Desktop/Single_cell_output/plot_integrated.svg", width = 3.2,   height = 2, dpi ="retina")

sc_nonint <- seur_obj %>%
  DimPlot(
    group.by = "seurat_clusters",
    label = T,
    label.size = 2.5,
    
  ) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(axis.text = element_blank()) +
  NoAxes()

x <- (nonint_sample +int_sample ) 

y <-  (nonint_injury+int_injury  )

fig1 <- (y|x)/ sc_int 
ggsave("/Users/ranojoychatterjee/Desktop/Single_cell_output/plot_integrated.png", width = 6.8,   height = 4, dpi ="300")


plots <- align_plots(p3, p1, align = 'v', axis = 'l')
# then build the bottom row
bottom_row <- plot_grid(plots[[2]], p2, labels = c('B', 'C'), label_size = 12)

# then combine with the top row for final plot
plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)


bottom_row <- ggplot(sc_int)

seur_obj <- FindVariableFeatures(seur_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

BuildClusterTree(object = seur_obj, 
                 genes.use =seur_obj@assays$integrated@var.features, 
                 do.plot = TRUE, 
                 do.reorder = TRUE, 
                 reorder.numeric = TRUE,
                 show.progress = TRUE)


all.markers <- FindAllMarkers(object = seur_obj)
top10 <- all.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = obj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)


DoHeatmap(subset(seur_obj, downsample = 100), features = features, size = 3)




seur_obj$Celltype <- sub('Fibroblasts 1', "KC", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 2", "SG1", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 3", "SG2", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 4", "Mac", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 5", "Mast", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 6", "FB 7", seur_obj$Celltype)
seur_obj$Celltype <- sub("Fibroblasts 7", "APM", seur_obj$Celltype)
seur_obj$Celltype <- sub("Keratinocytes", "Mono", seur_obj$Celltype)
seur_obj$Celltype <- sub("Adipocyte", "Mu 2", seur_obj$Celltype)
seur_obj$Celltype <- sub("Sebaceous gland cells 1", "FB 2", seur_obj$Celltype)
seur_obj$Celltype <- sub("Sebaceous gland cells 2", "EC 1", seur_obj$Celltype)
seur_obj$Celltype <- sub("SK Muscle 1", "Adipo", seur_obj$Celltype)
seur_obj$Celltype <- sub("SK Muscle 2", "Mu 1", seur_obj$Celltype)
seur_obj$Celltype <- sub("SK Muscle 3", "Mu 3", seur_obj$Celltype)
seur_obj$Celltype <- sub("Endothelial 1", "FB 3", seur_obj$Celltype)
seur_obj$Celltype <- sub("Endothelial 2", "FB 5", seur_obj$Celltype)
seur_obj$Celltype <- sub("Errector Pilli", "FB 4", seur_obj$Celltype)
seur_obj$Celltype <- sub("Immune", "EC 2", seur_obj$Celltype)
seur_obj$Celltype <- sub("Macrophages 1", "FB 1", seur_obj$Celltype)
seur_obj$Celltype <- sub("Macrophages 2", "FB 6", seur_obj$Celltype)
seur_obj$Celltype <- sub("Monocytes", "Macro", seur_obj$Celltype)
