### Monocle3 ###

library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)

setwd("C:/Users/ranoj/Desktop/Single_Cell_output/analysis")


sample <-"NonWound1"

cds <- load_cellranger_data(paste0("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/", sample, sep = ""))

cds <- preprocess_cds(cds, num_dim = 100)

plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds)

plot_cells(cds)

cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds)
#plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="cluster_row_col",
                    max.size=3)
assigned_type_marker_test_res <- top_markers(cds,
                                             group_cells_by="cluster",
                                             reference_cells=1000,
                                             cores=8)

garnett_markers <- assigned_type_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)
# Exclude genes that are good markers for more than one cell type:
garnett_markers <- garnett_markers %>% 
  group_by(gene_short_name) %>%
  filter(n() == 1)

monocle_object <- cds
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object,reduction_method = "UMAP")
non_wound<-plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE)+ theme(legend.position = "none")

cds_subset <- choose_cells(cds)

sample <-"Wound1"
cds_wound <- load_cellranger_data(paste0("C:/Users/ranoj/Box/snRNA_CellRanger_Wound_nonWound/data/", sample, sep = ""))

cds_wound <- preprocess_cds(cds_wound , num_dim = 100)

plot_pc_variance_explained(cds_wound )

cds_wound  <- reduce_dimension(cds_wound )

plot_cells(cds_wound )

cds_wound  <- cluster_cells(cds_wound , resolution=1e-5)
plot_cells(cds_wound )

plot_genes_by_group(cds_wound,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=3)




monocle_object <- cds_wound
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object,reduction_method = "UMAP")
wound<-plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE) + theme(legend.position = "none")
cds_wound_subset <- choose_cells(cds_wound)

cds_combine <- combine_cds(list(cds_subset, cds_wound_subset))

cds_combine <- preprocess_cds(cds_combine , num_dim = 100)
cds_combine <- reduce_dimension(cds_combine )
cds_combine  <- cluster_cells(cds_combine , resolution=1e-5)
plot_cells(cds_combine )
monocle_object <- cds_combine
monocle_object <- learn_graph(monocle_object, use_partition = TRUE)
monocle_object <- order_cells(monocle_object,reduction_method = "UMAP")
combined <- plot_cells(monocle_object,
           color_cells_by = "pseudotime",
           graph_label_size=5,
           show_trajectory_graph = TRUE)  + 
  theme(legend.title=element_text(size=6),
                                                  legend.text=element_text(size=5),
                                                  legend.position = c(0.85, 0.3))

fig_plot <- ggdraw() +
  draw_plot(non_wound_try,
            x = 0, width = .40) +
  draw_plot(wound_try + theme(axis.title.y = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.line.y = element_blank()),
            x = .40, width = .25) +
  draw_plot(combined + theme(axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.line.y = element_blank()),
            x = .65, width = .365) +
  draw_plot_label((label = c("a", "Normal", "Injured" ,"Normal VS Injured" )),
                  size = 7,x = c(0, .18, .45, .70), y = c(1, 0.99, 0.99,0.99 ))

ggsave(filename = 'trajectory.pdf', plot=fig_plot, width = 181, height = 110, units = "mm", dpi = 350)
ggsave(filename = 'trajectory.png', plot=fig_plot, width = 181, height = 110, units = "mm", dpi = 350)
