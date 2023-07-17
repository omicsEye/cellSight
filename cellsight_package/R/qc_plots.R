
#' Quality control filtering with plots
#'
#' @param seur_object The seurat object from the previous directory step
#'
#' @return QC plots
#' @export
#'
#' @examples

qc_plot <- function(seur_object){
  obj_list <- seur_object
  for(index in 1:length(obj_list)){
    print(index)
    obj_list[[index]][["percent.mt"]] <- obj_list[[index]] %>%
      PercentageFeatureSet(pattern = "^mt-")

    plot_violin <- VlnPlot(obj_list[[index]],features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                           ncol = 3)
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_violin_",index,".png"))
    plot_scatter <- FeatureScatter(obj_list[[index]],feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_scatter_",index,".png"),bg="white")
    grid <- plot_grid(
      plot_violin, plot_scatter,
      labels = "AUTO", ncol = 1
    )
    ggsave(paste0("/Users/ranoj/Desktop/Single_cell_output/qc/qcplot_grid_joined_",index,".png"),bg="white")
  }
}
