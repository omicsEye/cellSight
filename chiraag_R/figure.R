library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr) 
library(omicsArt)
library("ggplot2")
library(cowplot)

seur_obj <-
  readRDS("C:/Users/ranoj/Desktop/Single_cell_output/objects/final_integrated_clustered1118.rds")

new.cluster.ids <- c("Fibroblasts 1", "Fibroblasts 2", "Fibroblasts 3", "Keratinocytes", "Fibroblasts 4", "Monocytes",
                     "Fibroblasts 5", "Macrophages 1", "Sebaceous gland cells 1","Endothelial","Errector Pilli","Endothelial","Macrophages 2",
                     "Fibroblasts 6","Sk Mucle 1","SK Muscle 2","Adipocyte","Sk Muscle 3","Fibroblasts 7","Sebaceous gland cells 1", "Immune")
names(new.cluster.ids) <- levels(seur_obj)
seur_obj$Celltypes <- RenameIdents(seur_obj, new.cluster.ids)


cells.use <- WhichCells(seur_obj, idents = 'NA')
seur_obj <- SetIdent(seur_obj, cells = cells.use, value = 'Sebaceous gland cells 1')

cells.use <- WhichCells(seur_obj, idents = 'Unknown 1')
seur_obj <- SetIdent(seur_obj, cells = cells.use, value = 'Sebaceous gland cells 2')

cells.use <- WhichCells(seur_obj, idents = 'Unknown 2')
seur_obj <- SetIdent(seur_obj, cells = cells.use, value = 'Immune')


seur_obj$Celltype <- Idents(seur_obj)
seur_obj$type <- seur_obj$type %>%
  as.factor()

lapply(
  levels(seur_obj[["Celltype"]][[1]]),
  function(x)FindMarkers(seur_obj,ident.1 = x,min.pct = 0.25)
) -> cluster.markers


# This simply adds the cluster number to the results of FindMarkers
sapply(0:(length(cluster.markers)-1),function(x) {
  cluster.markers[[x+1]]$gene <<- rownames(cluster.markers[[x+1]])
  cluster.markers[[x+1]]$cluster <<- x
})

cluster.markers <- as_tibble(do.call(rbind,cluster.markers)) %>% arrange(p_val_adj) 
cluster.markers

cluster.markers %>%
  group_by(cluster) %>%
  slice(5) %>%
  pull(gene) -> best.wilcox.gene.per.cluster

best.wilcox.gene.per.cluster

VlnPlot(seur_obj,features=best.wilcox.gene.per.cluster)

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = 7, angle = 0), 
          axis.text.y = element_text(size =5), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(size =7,angle = 45), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


#imp_genes <- c( "Adgre1", "Ptprc","Pdgfra", "Pdgfrb", "Col1a1", "Pparg","Krt77","Pecam1","Cdh4","Krt77","Cd44","Adipor2","Ppard","Cd36","Cd55","Cd46","Fap")

#imp_genes <- c( "Gm26917","Tgfbi","Gm42418","Skint6","Ebf1","Ciita","Gm42418","Nrros","Auts2","Tshz2","Myh11","Mecom","Spag17","Ttn","Slc1a3","Cmss1","Gpc5","Far2","St6galnac3")

imp_genes <- c("Col1a1","Fbln2","Adgre1","Ptprc","Itga8","Tnnt3","Pcdh9","Flt1","Ebf1","Nrros","Auts2","Tshz2","Myh11","Mecom","Spag17","Ttn","Slc1a3","Cmss1","Gpc5","Far2","St6galnac3")

p <- StackedVlnPlot(obj = seur_obj, features = imp_genes)

p
png(paste0(fig_dir, "new_feature_plot.png"), width=7.2, height=2.75, res=250, units='in')
p 
dev.off()
#p + theme(axis.text.y = element_text(face = "bold", color = "blue",size = 12, angle = 45))

png(paste0(fig_dir, "new_feature_plot.png"), width=7, height=5.5, res=250, units='in')
#(p1|p2)/p3 +umap_theme
#CombinePlots(, ncol=3)
p
dev.off()


library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(viridis)

umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

fig_dir <- "C:/Users/ranoj/Desktop/Single_cell_output/figure/"

pdf(paste0(fig_dir, "umap_clusters_sct.pdf"), width=2.4, height=2.75)
p1 <- DimPlot(seur_obj, group.by='type', reduction='umap',raster = 10000) +
  ggtitle('Type') + umap_theme +  theme(plot.title = element_text(size = 7, face = "bold"))+
  scale_color_manual(labels = c("Unwounded", "Wounded"), values = c("green", "red")) +
  theme(legend.text=element_text(size=5))


print(p1)

dev.off()

pdf(paste0(fig_dir, "umap_batch_sct.pdf"), width=2.4, height=2.75)
p2 <- DimPlot(seur_obj, group.by='sample', reduction='umap',raster = 10000) +
  ggtitle('Sample') + umap_theme+  theme(plot.title = element_text(size = 7, face = "bold"))+
  scale_color_manual(labels = c("Unwounded_1","Unwounded_2", "Wounded_1","Wounded_2"), values = c("#0CB702","#00C19A" ,"#FF68A1","#F8766D")) +
  theme(legend.text=element_text(size=5))
print(p2)
dev.off()

pdf(paste0(fig_dir, "umap_cluster_sct.png"), width=5, height=2.75, res=250, units='in')
p3 <- DimPlot(seur_obj, group.by='Celltype', reduction='umap',label = T, label.size = 3,repel = T,raster = 10000) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none") +theme_nothing() +
  theme(axis.text = element_blank())+
  ggtitle('Seurat Clusters') + umap_theme +  theme(plot.title = element_text(size = 7, face = "bold"))+
  theme(legend.text=element_text(size=5))

png(paste0(fig_dir, "umap_cluster_sct.png"), width=5, height=2.75, res=250, units='in')

p3
dev.off()

png(paste0(fig_dir, "pngs/umap_joined_all.png"), width=9.8, height=2.75, res=250, units='in')
(p1|p2|p3) +umap_theme
#CombinePlots(, ncol=3)
dev.off()

box_association <- readRDS("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/analysis/tweedieverse_new/cluster_Fibroblasts 1/figures/type_gg_associations.rds")

fig2 <- ggdraw() + draw_plot(box_association[[22]] )
            